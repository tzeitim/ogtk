import os
import json
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Set, Iterable
import polars as pl

from textual.app import App, ComposeResult
from textual.widgets import DirectoryTree, Tree, Footer, Header, Static, Button, Label, SelectionList, Log
from textual.containers import Horizontal, Vertical, Container
from textual.screen import Screen, ModalScreen
from textual.reactive import var
from textual.binding import Binding
from textual.message import Message

from rich.json import JSON
from rich.panel import Panel
from rich.table import Table


from ogtk.ltr.fracture.post.metrics.summary import PipelineMetricsCollection

import re
from datetime import datetime
import logging


pl.Config.set_float_precision(2)
pl.Config.set_thousands_separator("_")

def sanitize_id(name: str) -> str:
    sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', name.lower())
    if sanitized and sanitized[0].isdigit():
        sanitized = f"id_{sanitized}"
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_')
    return sanitized

def _format_number(value, metric_name) -> str:
    """Format numeric values with underscore separators for thousands."""
    if isinstance(value, (int, float)):
        # For percentages, rates, and similar metrics, format with 2 decimal places
        if isinstance(value, float):
            if any(keyword in metric_name.lower() for keyword in ["fraction", "rate", "mean", "median", "percent"]):
                return f"{value:.2f}"
            elif abs(value) >= 1000:
                return f"{value:,.2f}".replace(",", "_")
            else:
                return f"{value:.2f}"  
        
        # For large integers, add underscore separators
        if isinstance(value, int) and abs(value) >= 1000:
            return f"{value:_}".replace(",", "_")
        
        # For large floats (non-percentages), add underscore separators
        if isinstance(value, float) and abs(value) >= 1000:
            int_part = int(value)
            decimal_part = value - int_part
            formatted_int = f"{int_part:_}".replace(",", "_")
            return f"{formatted_int}{str(decimal_part)[1:]}"
            
    # Return as string for non-numeric or small numbers
    return str(value)

class SampleSelected(Message):
    """Message sent when a sample is selected."""
    
    def __init__(self, sample_path: Path):
        super().__init__()
        self.sample_path = sample_path


class FilteredDirectoryTree(DirectoryTree):
    def filter_paths(self, paths: Iterable[Path]) -> Iterable[Path]:
        """Filter to only show directories, JSON files, and Markdown files."""
        # This must return an iterable, not a list comprehension directly
        for path in paths:
            # Always include directories
            if path.is_dir() and not path.name.startswith(".") and (path/'pipeline_summary.json').exists() :
                yield path
            # Only include .json and .md files
            elif path.suffix.lower() in (".json", ".md"):
                yield path


class JSONDisplay(Static):
    """Widget to display JSON data."""
    
    def __init__(
        self,
        json_data: Optional[Dict[str, Any]] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        classes: Optional[str] = None,
    ):
        super().__init__(name=name, id=id, classes=classes)
        self.json_data = json_data or {}
    
    def update_data(self, json_data: Dict[str, Any]) -> None:
        """Update the JSON data."""
        self.json_data = json_data
        self.update(self._render_json())
    
    def _render_json(self) -> str:
        """Render the JSON data."""
        if not self.json_data:
            return "No data selected"
        
        return JSON.from_data(self.json_data)
    
    def compose(self) -> ComposeResult:
        self.update(self._render_json())
        yield from []


class FigureViewer(ModalScreen):
    """Modal screen to view figure files."""

    BINDINGS = [
        Binding("escape", "close", "Close"),
        Binding("q", "close", "Close"),
    ]

    CSS = """
    #figure-container {
        width: 100%;
        height: 100%;
        align: center middle;
    }

    Static {
        width: 100%;
        height: auto;
        content-align: center middle;
    }

    .figure-button {
        margin: 1;
    }
    """

    def __init__(self, figure_paths: List[str], experiment_dir: Path):
        """Initialize with figure paths."""
        super().__init__()
        self.figure_paths = figure_paths
        self.experiment_dir = experiment_dir
        self.current_figure = None

    def action_close(self) -> None:
        """Close the modal."""
        self.app.pop_screen()

    def show_figure(self, path: str) -> None:
        """Show the selected figure."""
        self.current_figure = path
        self.notify(f"Viewing {os.path.basename(path)}")

        # Try to open the image in an external viewer
        try:
            abs_path = self._resolve_figure_path(path)
            if sys.platform == "darwin":  # macOS
                subprocess.Popen(["open", abs_path])
            elif sys.platform == "win32":  # Windows
                os.startfile(abs_path)
            else:  # Linux
                subprocess.Popen(["xdg-open", abs_path])
        except Exception as e:
            self.notify(f"Error opening figure: {str(e)}", severity="error")

    def _resolve_figure_path(self, path: str) -> str:
        """Resolve the figure path relative to the experiment directory."""
        if path.startswith("../"):
            # Handle relative paths in the JSON
            return str(self.experiment_dir.parent / path[3:])
        return str(self.experiment_dir / path)

    def compose(self) -> ComposeResult:
        """Compose the screen layout."""
        yield Header(show_clock=True)
        with Container(id="figure-container"):
            with Vertical():
                with Horizontal():
                    for i, path in enumerate(self.figure_paths):
                        filename = os.path.basename(path)
                        yield Button(filename, classes="figure-button", id=f"fig-{i}")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button press event."""
        button_id = event.button.id
        if button_id and button_id.startswith("fig-"):
            try:
                idx = int(button_id.split("-")[1])
                if 0 <= idx < len(self.figure_paths):
                    self.show_figure(self.figure_paths[idx])
            except (ValueError, IndexError):
                self.notify("Invalid figure index", severity="error")
###
class SmartExperimentDirectoryTree(DirectoryTree):
    """Directory tree that adapts to workdir or single experiment structure."""
    BINDINGS = [
        Binding("j", "cursor_down", "Down"),
        Binding("k", "cursor_up", "Up"),
        Binding("enter", "select_cursor", "Select"),
        Binding("space", "select_cursor", "Select"),
    ]
    
    def __init__(self, path: Path, **kwargs):
        super().__init__(path, **kwargs)
        self.selected_samples = set()
        self.sample_nodes = {}  # Track sample nodes for styling

    
    def filter_paths(self, paths: Iterable[Path]) -> Iterable[Path]:
        """Filter to show experiment directories and sample directories."""
        for path in paths:
            if not path.is_dir() or path.name.startswith('.'):
                continue
                
            # Show sample directories with pipeline_summary.json
            if self._is_valid_sample_dir(path):
                yield path
            # Show experiment directories that contain valid samples
            elif self._is_experiment_dir(path):
                yield path

    def deselect_all_samples(self) -> None:
       """Deselect all selected sample directories."""
       for sample_id in list(self.selected_samples):
           if sample_id in self.sample_nodes:
               self._update_node_style(self.sample_nodes[sample_id], selected=False)
       
       self.selected_samples.clear()
       self.sample_nodes.clear()

    def select_all_samples(self) -> int:
        """Select all visible sample directories."""
        selected_count = 0
        
        def visit_node(node):
            nonlocal selected_count
            if hasattr(node, 'data') and hasattr(node.data, 'path'):
                if self._is_valid_sample_dir(node.data.path):
                    sample_id = self._get_sample_id(node.data.path)
                    if sample_id not in self.selected_samples:
                        self.selected_samples.add(sample_id)
                        self._update_node_style(node, selected=True)
                        self.sample_nodes[sample_id] = node
                        selected_count += 1
                        # this makes the selection visible by the main app
                        self.post_message(SampleSelected(node.data.path))
            
            for child in node.children:
                visit_node(child)
        
        visit_node(self.root)
        return selected_count

    def _is_experiment_dir(self, path: Path) -> bool:
        """Check if this directory contains sample subdirectories."""
        try:
            # Look for subdirectories that contain pipeline_summary.json
            for subdir in path.iterdir():
                if subdir.is_dir() and (subdir / "pipeline_summary.json").exists():
                    return True
        except PermissionError:
            pass
        return False
    
    def _is_valid_sample_dir(self, path: Path) -> bool:
        """Check if this is a valid sample directory."""
        return (path / "pipeline_summary.json").exists()

    def on_tree_node_highlighted(self, event: Tree.NodeHighlighted) -> None:
        """Handle node highlighting for live display updates."""
        if not hasattr(event.node, 'data') or not hasattr(event.node.data, 'path'):
            return
 
        selected_path = event.node.data.path
        # Only process sample directories for live display
        if self._is_valid_sample_dir(selected_path):
            # Post message for live display update
            self.post_message(SampleSelected(selected_path))

    def on_tree_node_selected(self, event: Tree.NodeSelected) -> None:
        """Handle node selection - works with Textual's Tree widget."""
        if not hasattr(event.node, 'data') or not hasattr(event.node.data, 'path'):
            return
            
        selected_path = event.node.data.path
        
        # Only process sample directories
        if self._is_valid_sample_dir(selected_path):
            sample_id = self._get_sample_id(selected_path)
            
            # Toggle selection
            if sample_id in self.selected_samples:
                self.selected_samples.remove(sample_id)
                self._update_node_style(event.node, selected=False)
                #self.app.notify(f"Sample '{sample_id}' deselected")
            else:
                self.selected_samples.add(sample_id)
                self._update_node_style(event.node, selected=True)
                #self.app.notify(f"Sample '{sample_id}' selected")
            
            # Store node reference for future styling updates
            self.sample_nodes[sample_id] = event.node
            
            # Trigger sample selection event for the main app
            self.post_message(SampleSelected(selected_path))
    
    def _update_node_style(self, node, selected: bool) -> None:
        """Update the visual style of a node based on selection state."""
        if selected:
            #node.add_class("selected-sample")
            # Update the label to show selection
            if hasattr(node, '_label'):
                original_label = str(node._label).replace(" ✓", "")
                node.set_label(f"{original_label} ✓")
        else:
            #node.remove_class("selected-sample")
            # Remove selection indicator from label
            if hasattr(node, '_label'):
                original_label = str(node._label).replace(" ✓", "")
                node.set_label(original_label)
    
    def _get_sample_id(self, sample_path: Path) -> str:
        """Get the sample ID based on structure type."""
        return f"{sample_path.parent.name}/{sample_path.name}"
    
    def get_selected_sample_paths(self) -> list[Path]:
        """Get paths to all selected samples."""
        selected_paths = []
        for sample_id in self.selected_samples:
            # Sample ID is always in experiment/sample format
            experiment_name, sample_name = sample_id.split('/', 1)
            sample_path = self.path / experiment_name / sample_name
            
            if sample_path.exists():
                selected_paths.append(sample_path)
        return selected_paths
    
    def get_selected_sample_ids(self) -> set[str]:
        """Get the selected sample IDs."""
        return self.selected_samples.copy()
    
    def clear_selection(self) -> None:
        """Clear all selected samples."""
        for sample_id in list(self.selected_samples):
            if sample_id in self.sample_nodes:
                self._update_node_style(self.sample_nodes[sample_id], selected=False)
        
        self.selected_samples.clear()
        self.sample_nodes.clear()
        self.app.notify("All samples deselected")
    
    def select_sample_by_id(self, sample_id: str) -> bool:
        """Programmatically select a sample by ID. Returns True if successful."""
        if '/' not in sample_id:
            return False
        
        # Sample ID is always in experiment/sample format
        experiment_name, sample_name = sample_id.split('/', 1)
        sample_path = self.path / experiment_name / sample_name
        
        if not sample_path.exists() or not self._is_valid_sample_dir(sample_path):
            return False
        
        if sample_id not in self.selected_samples:
            self.selected_samples.add(sample_id)
            
            # Find and update the corresponding node if it exists
            for node_id, node in self.nodes.items():
                if (hasattr(node, 'data') and 
                    hasattr(node.data, 'path') and 
                    node.data.path == sample_path):
                    self._update_node_style(node, selected=True)
                    self.sample_nodes[sample_id] = node
                    break
            
            return True
        
        return False
###
class ExperimentList(Static):
    """Widget that displays a list of valid experiments"""
    def __init__(
        self,
        path: Path,
        id: Optional[str] = None,
        name: Optional[str] = None,
        classes: Optional[str] = None,
    ):
        super().__init__(name=name, id=id, classes=classes)
        self.path = path
   
    def _populate_list(self) -> list[tuple[str, str]]:
        """Filter to only show directories, JSON files, and Markdown files."""
        files = []
        for path in self.path.iterdir():
            if path.is_dir() and not path.name.startswith(".") and (path/'pipeline_summary.json').exists() :
                files.append(path)
        return [(f.name, f.name) for f in files]

    def compose(self) -> ComposeResult:
        yield SelectionList[str](*self._populate_list())


class MetricsTable(Static):
    """Widget to display metrics as a table."""

    def __init__(
        self,
        metrics: Optional[Dict[str, Any]] = None,
        stage: str = "",
        name: Optional[str] = None,
        id: Optional[str] = None,
        classes: Optional[str] = None,
    ):
        super().__init__(name=name, id=id, classes=classes)
        self.metrics = metrics or {}
        self.stage = stage

    def update_metrics(self, metrics: Dict[str, Any], stage: str, path: Path|None = None) -> None:
        """Update the metrics data."""
        self.metrics = metrics
        self.stage = stage
        self.path = path
        self.update(self._render_metrics())

    def _format_number(self, value, metric_name) -> str:
        """Format numeric values with underscore separators for thousands."""
        return _format_number(value, metric_name)

    def _render_metrics(self) -> Panel:
        """Render the metrics as a table."""
        if not self.metrics:
            return Panel("No metrics selected")

        table = Table(title=self.path.name)
        table.add_column("Metric")
        table.add_column("Value")

        for key, value in self.metrics.items():
            formatted_value = self._format_number(value, key)
            table.add_row(key, formatted_value)

        return Panel(table, title=f"{self.path.as_posix()} Metrics")


    def compose(self) -> ComposeResult:
        self.update(self._render_metrics())
        yield from []


class ComparisonDataTable(Static):
    """Widget to display comparison data in a table."""

    def __init__(
        self,
        df: Optional[pl.DataFrame] = None,
        title: str = "Sample Comparison",
        name: Optional[str] = None,
        id: Optional[str] = None,
        classes: Optional[str] = None,
    ):
        super().__init__(name=name, id=id, classes=classes)
        self.df = df or pl.DataFrame()
        self.title = title

    def update_data(self, df: pl.DataFrame, title: str = None) -> None:
        """Update the dataframe data."""
        self.df = df
        if title:
            self.title = title
        self.update(self._render_table())

    def _render_table(self) -> Panel:
        """Render the dataframe as a table."""
        if self.df.is_empty():
            return Panel("No comparison data available")

        table = Table(title=self.title)

        # Add columns
        for column in self.df.columns:
            table.add_column(column)

        # Add rows
        for row in self.df.to_dicts():
            formatted_row = []
            for col in self.df.columns:
                value = row[col]
                if col == "sample_id":  
                    formatted_row.append(str(value))
                else:
                    formatted_row.append(_format_number(value, col))
            table.add_row(*formatted_row)

        return Panel(table, title=self.title)

    def compose(self) -> ComposeResult:
        self.update(self._render_table())
        yield from []


class ComparisonScreen(Screen):
    """Screen for comparing metrics across samples."""

    CSS_PATH = "comparison_screen.tcss"

    BINDINGS = [
        Binding("q", "return_to_main", "Return to Main"),
        Binding("r", "run_comparison", "Run Comparison"),
    ]

    def __init__(self, collection: PipelineMetricsCollection, preselected_samples=None):
        super().__init__()
        self.collection = collection
        self.selected_samples = set(preselected_samples or [])

    def action_return_to_main(self) -> None:
        """Return to main screen."""
        self.app.pop_screen()

    def action_run_comparison(self) -> None:
        """Run the comparison with selected samples."""
        self._run_comparison()

    def action_focus_next_button(self) -> None:
        """Focus the next button, skipping checkboxes."""
        buttons = list(self.query("Button"))
        if not buttons:
            return
        
        current_focus = self.focused
        if current_focus and current_focus in buttons:
            current_index = buttons.index(current_focus)
            next_index = (current_index + 1) % len(buttons)
            buttons[next_index].focus()
        else:
            buttons[0].focus()

    def action_focus_previous_button(self) -> None:
        """Focus the previous button, skipping checkboxes."""
        buttons = list(self.query("Button"))
        if not buttons:
            return
        
        current_focus = self.focused
        if current_focus and current_focus in buttons:
            current_index = buttons.index(current_focus)
            prev_index = (current_index - 1) % len(buttons)
            buttons[prev_index].focus()
        else:
            # If not currently on a button, focus the last button
            buttons[-1].focus()

    def compose(self) -> ComposeResult:
        """Compose the screen layout."""
        yield Header(show_clock=True)

        with Vertical(id="left-panel"):
            with Container(id="sample-selection"):
                yield Label("Pre-selected samples to compare:", id="selection-label")

                selection_options = [(sample_id, sample_id) for sample_id in self.selected_samples]
                selection_list = SelectionList[str](*selection_options, id="sample-selection-list")
                
                # Pre-select all items
                for sample_id in self.selected_samples:
                    selection_list.select(sample_id)
                
                yield selection_list

            with Vertical(id="metric-buttons"):
                yield Button("Compare Valid UMI %", classes="metric-button", id="compare-valid-umi")
                yield Button("Compare Read Coverage", classes="metric-button", id="compare-read-coverage")

        with Vertical(id="right-panel"):
            yield ComparisonDataTable(id="comparison-table")

        yield Footer()

    def on_mount(self) -> None:
        """ """
        self.call_after_refresh(self._run_comparison) 

    def _run_comparison(self) -> None:
        """Run the selected comparison."""
        selection_list = self.query_one("#sample-selection-list", SelectionList)
        selected_indices = selection_list.selected
        
        if not selected_indices:
            self.notify("Please select at least one sample to compare", severity="warning")
            return

        current_selected = list(self.selected_samples)

        df = self.collection.as_wide_df(strip_step_names=True)
        df = df.filter(pl.col("sample_id").is_in(list(current_selected)))

        if df.height == 0:
            self.notify("No valid data for selected samples", severity="warning")
            return

        self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Valid UMI Comparison")


    def _get_selected_sample_ids(self) -> Set[str]:
        """Get currently selected sample IDs from SelectionList."""
        selection_list = self.query_one("#sample-selection-list", SelectionList)
        return set(selection_list.selected)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        button_id = event.button.id

        if button_id == "compare-button":
            self._run_comparison()

        elif button_id == "compare-valid-umi":
            current_selected = self._get_selected_sample_ids()

            if not current_selected:
                self.notify("Please select at least one sample to compare", severity="warning")
                return

            df = self.collection.get_valid_umi_stats()
            df = df.filter(pl.col("sample_id").is_in(list(current_selected)))

            if df.height == 0:
                self.notify("No valid UMI data for selected samples", severity="warning")
                return

            self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Valid UMI Comparison")

        elif button_id == "compare-read-coverage":
            current_selected = self._get_selected_sample_ids()

            if not current_selected:
                self.notify("Please select at least one sample to compare", severity="warning")
                return

            df = self.collection.calculate_read_coverage()
            df = df.filter(pl.col("sample_id").is_in(list(current_selected)))

            if df.height == 0:
                self.notify("No read coverage data for selected samples", severity="warning")
                return

            self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Read Coverage Comparison")

                

class FractureExplorer(App):
    """A Textual app to explore Fracture pipeline outputs."""
    TITLE = """ ○━━━━━━━━○ """
    OK = """
    ███████╗██████╗  █████╗  ██████╗████████╗██╗   ██╗██████╗ ███████╗
    ██╔════╝██╔══██╗██╔══██╗██╔════╝╚══██╔══╝██║   ██║██╔══██╗██╔════╝
    █████╗  ██████╔╝███████║██║        ██║   ██║   ██║██████╔╝█████╗  
    ██╔══╝  ██╔══██╗██╔══██║██║        ██║   ██║   ██║██╔══██╗██╔══╝  
    ██║     ██║  ██║██║  ██║╚██████╗   ██║   ╚██████╔╝██║  ██║███████╗
    ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝
    """

    CSS_PATH = "fracture_viewer.tcss"

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("r", "refresh", "Refresh"),
        Binding("f", "view_figures", "View Figures"),
        Binding("c", "compare_samples", "Compare Samples"),
        Binding("e", "toggle_files", "Toggle Files"),
        Binding("a", "select_all", "Select All Files"),
        Binding("l", "toggle_log", "Toggle Log"),  
        Binding("d", "deselect_all", "Clear Selection"),
    ]
    show_tree = var(False)
    show_log = var(False)


    def action_deselect_all(self) -> None:
        """Select All Files"""
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        tree.deselect_all_samples()

    def action_select_all(self) -> None:
        """Select All Files"""
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        selected_n = tree.select_all_samples()
        self.notify(f"Selected {selected_n} samples")

    def action_toggle_log(self) -> None:
        """Called in response to log toggle key binding."""
        self.show_log = not self.show_log

    def watch_show_log(self, show_log: bool) -> None:
        """Called when show_log is modified."""
        self.set_class(show_log, "-show-log")

    def action_toggle_files(self) -> None:
        """Called in response to key binding."""
        self.show_tree = not self.show_tree
    
    def watch_show_tree(self, show_tree: bool) -> None:
        """Called when show_tree is modified."""
        self.set_class(show_tree, "-show-raw")

    def log_message(self, message: str, level: str = "INFO") -> None:
        """Write a message to the log widget."""
        log = self.query_one("#log-display", Log)
        timestamp = datetime.now().strftime("%H:%M:%S")
        log.write_line(f"[{timestamp}] {level}: {message}")

    def __init__(self, experiment_dir: str):
        """Initialize the app with the experiment directory."""
        super().__init__()
        # Initialize log buffer for early messages
        self.log_buffer = []
        self.log_widget_ready = False

        # Log the very first message
        timestamp = datetime.now().strftime("%H:%M:%S")
        first_msg = f"[{timestamp}] INFO: App initialization started"
        self.log_buffer.append(first_msg)
        print(first_msg)  # Also print to terminal immediately
        
        # Setup early logging to buffer
        self._setup_early_logging()
        
        self.experiment_dir = Path(experiment_dir).resolve()
        if not self.experiment_dir.exists() or not self.experiment_dir.is_dir():
            raise ValueError(f"Experiment directory {experiment_dir} does not exist.")

        self.current_figures = []
        self.current_sample_dir = None

        # For sample comparison
        self.selected_samples = set()
        self.metrics_collection = PipelineMetricsCollection()
        self.pipeline_files = []

    def _setup_early_logging(self):
        """Setup logging to buffer early messages."""
        # Configure logging to go to both buffer and stdout
        class BufferHandler(logging.Handler):
            def __init__(self, app):
                super().__init__()
                self.app = app
                
            def emit(self, record):
                msg = self.format(record)
                self.app.log_buffer.append(msg)
        
        # Setup dual logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s [%(levelname)s] %(message)s',
            handlers=[
                logging.StreamHandler(sys.stdout),
                BufferHandler(self)
            ]
        )
    def hybrid_log(self, message: str, level: str = "INFO") -> None:
        """Log to buffer early, then to widget once available."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        formatted_msg = f"[{timestamp}] {level}: {message}"
        
        if self.log_widget_ready:
            # Widget is ready, log directly
            try:
                log = self.query_one("#log-display", Log)
                log.write_line(formatted_msg)
            except:
                # Fallback to buffer if widget query fails
                self.log_buffer.append(formatted_msg)
        else:
            # Buffer the message
            self.log_buffer.append(formatted_msg)
            # Also log to terminal
            #logging.info(message)

    def compose(self) -> ComposeResult:
        """Compose the app layout."""
        self.hybrid_log("Starting UI composition", "INFO")

        
        #yield Static(title_art, id="title-art")
        yield Header(show_clock=True)
        
        with Vertical(id='experiments'):
            self.hybrid_log("Initializing directory tree scanner...", "INFO")
            yield SmartExperimentDirectoryTree(self.experiment_dir, id="experiment-tree")

        with Vertical(id="metrics-container"):
            yield MetricsTable(id="all-metrics", stage="all")
        
            # Add a JSONDisplay to show the full JSON data
            yield JSONDisplay(id="raw-json-display")
            yield Log(id="log-display")

            with Horizontal(id="control-buttons"):
                yield Button("Select All", id="select-all-button")
                yield Button("Clear Selection", id="deselect-all-button")
                yield Button("Compare Samples", id="compare-samples-button")
                yield Button("Select Sample", id="select-sample-button")
                yield Button("View Figures", id="figures-button")

        yield Footer()

    def on_mount(self) -> None:
        """Setup when the app is mounted."""
        self.theme = "gruvbox"
        
        self._initialize_log_widget()

        # Log when tree has finished loading
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        sample_count = len(list(tree.path.rglob("pipeline_summary.json")))
        self.hybrid_log(f"Directory tree loaded - found {sample_count} samples", "SUCCESS")

    def _initialize_log_widget(self) -> None:
        """Initialize log widget and replay buffered messages."""
        try:
            log = self.query_one("#log-display", Log)
            
            # Replay all buffered messages
            self.hybrid_log("=== LOG WIDGET INITIALIZED ===", "INFO")
            self.hybrid_log("Replaying buffered messages:", "INFO")
            
            for msg in self.log_buffer:
                log.write_line(msg)
            
            # Mark widget as ready
            self.log_widget_ready = True
            
            # Clear the buffer to save memory
            self.log_buffer.clear()
            
            self.hybrid_log("Log widget ready - future messages will appear here", "SUCCESS")
            
        except Exception as e:
            logging.error(f"Failed to initialize log widget: {e}")

    def on_sample_selected(self, event: SampleSelected) -> None:
        """Handle sample selection from the directory tree."""
        sample_dir = event.sample_path
        
        try:
            with open(sample_dir / "pipeline_summary.json", "r") as f:
                data = json.load(f)
            
            # Update the current sample directory
            self.current_sample_dir = sample_dir
            
            # Update select button text
            tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
            sample_id = tree._get_sample_id(sample_dir)
            select_button = self.query_one("#select-sample-button", Button)
            
            if sample_id in tree.selected_samples:
                select_button.label = f"Unselect '{sample_id}'"
            else:
                select_button.label = f"Select '{sample_id}'"
            
            # Reset figures list
            self.current_figures = []
            figures_button = self.query_one("#figures-button", Button)
            figures_button.remove_class("has-figures")
            
            # Update metrics tables
            is_parquet = "parquet" in data and "metrics" in data["parquet"]
            is_preprocess = "preprocess" in data and "metrics" in data["preprocess"]
            is_fracture = "fracture" in data and "metrics" in data["fracture"]
            is_all = is_parquet and is_preprocess and is_fracture

            if is_all:
                all_metrics = {}
                all_metrics.update(data["parquet"]["metrics"])
                all_metrics.update(data["preprocess"]["metrics"])
                all_metrics.update(data["fracture"]["metrics"])

                if "figures" in data["preprocess"]:
                    self.current_figures.extend(data["preprocess"]["figures"])
                    figures_button.add_class("has-figures")

                if "figures" in data["fracture"]:
                    self.current_figures.extend(data["fracture"]["figures"])
                    figures_button.add_class("has-figures")

                # this renders the JSON as table
                self.query_one("#all-metrics", MetricsTable).update_metrics(
                    all_metrics, stage="all", path=sample_dir
                )
            
            # Update JSON display
            self.query_one("#raw-json-display", JSONDisplay).update_data(data)
            
            # Check pipeline files for comparison
            summary_path = sample_dir / "pipeline_summary.json"
            if summary_path not in self.pipeline_files:
                self.pipeline_files.append(summary_path)
            
        except Exception as e:
            self.notify(f"Error loading JSON: {str(e)}", severity="error")

    def action_refresh(self) -> None:
        """Refresh the directory tree."""
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        tree.reload()

    def action_view_figures(self) -> None:
        """View figures for the current sample."""
        if self.current_figures and self.current_sample_dir:
            self.push_screen(FigureViewer(self.current_figures, self.current_sample_dir))
        else:
            self.notify("No figures available for the current sample", severity="warning")

    def action_compare_samples(self) -> None:
        """Open the sample comparison screen."""
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        selected_sample_ids = tree.get_selected_sample_ids()
        
        if not selected_sample_ids:
            self.notify("Please select at least one sample to compare", severity="warning")
            return
        
        # Update metrics collection from pipeline files
        if not self.pipeline_files:
            self.notify("No pipeline summary files loaded yet", severity="warning")
            return

        try:
            self.metrics_collection = PipelineMetricsCollection.from_summary_files(
                file_paths=self.pipeline_files
            )
            if not self.metrics_collection.samples:
                self.notify("No sample metrics found", severity="warning")
                return

            self.push_screen(ComparisonScreen(self.metrics_collection, preselected_samples=selected_sample_ids))
        except Exception as e:
            self.notify(f"Error loading metrics for comparison: {str(e)}", severity="error")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button press events."""
        button_id = event.button.id

        if button_id == "figures-button":
            self.action_view_figures()

        elif button_id == "select-all-button":
            self.action_select_all()

        elif button_id == "deselect-all-button":
            self.action_deselect_all()

        elif button_id == "compare-samples-button":
            self.action_compare_samples()

        elif button_id == "select-sample-button":
            if self.current_sample_dir:
                tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
                sample_id = tree._get_sample_id(self.current_sample_dir)
                select_button = self.query_one("#select-sample-button", Button)

                if sample_id in tree.selected_samples:
                    # Unselect the sample
                    tree.selected_samples.remove(sample_id)
                    if sample_id in tree.sample_nodes:
                        tree._update_node_style(tree.sample_nodes[sample_id], selected=False)
                        del tree.sample_nodes[sample_id]
                    
                    self.notify(f"Sample '{sample_id}' removed from selection")
                    select_button.label = f"Select '{sample_id}'"
                else:
                    # Select the sample
                    if tree.select_sample_by_id(sample_id):
                        self.notify(f"Sample '{sample_id}' added to selection")
                        select_button.label = f"Unselect '{sample_id}'"
                    else:
                        self.notify("Failed to select sample", severity="error")
            else:
                self.notify("No sample currently selected", severity="warning")

def main():
    """Console entry point for fracture-viewer."""
    import sys
    import os
    
    if len(sys.argv) > 1:
        experiment_dir = sys.argv[1]
    else:
        experiment_dir = os.environ.get("FRACTURE_EXPERIMENT_DIR", 
                                        os.path.expanduser("~/projects/lt/workdir"))
    
    app = FractureExplorer(experiment_dir)
    app.run()

if __name__ == "__main__":
    main()

#def main():
#    """Run the app."""
#    if len(sys.argv) > 1:
#        experiment_dir = sys.argv[1]
#    else:
#        experiment_dir = os.environ.get("FRACTURE_EXPERIMENT_DIR", 
#                                        os.path.expanduser("~/src/fracture-app/local-toy"))
#    
#    app = FractureExplorer(experiment_dir)
#    app.run()
#
#
#if __name__ == "__main__":
# 
#def main():
#    """Console entry point for fracture-viewer."""
#    import sys
#    import os
#    
#    if len(sys.argv) > 1:
#        experiment_dir = sys.argv[1]
#    else:
#        experiment_dir = os.environ.get("FRACTURE_EXPERIMENT_DIR", ".")
#    
#    # Import the app class (should be in the same file or imported)
#    app = FractureExplorer(experiment_dir)
#    app.run()
#
## Your existing FractureExplorer class here...
#
#if __name__ == "__main__":
#    main()
