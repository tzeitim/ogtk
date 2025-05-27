#!/usr/bin/env python
import os
import json
import sys
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Set, Iterable
import polars as pl
from datetime import datetime

from textual import on
from textual.app import App, ComposeResult
from textual.widgets import DirectoryTree, Tree, Footer, Header, Static, DataTable, Button, Markdown, Checkbox, Label, SelectionList, Pretty
from textual.containers import Horizontal, Vertical, Container, Grid
from textual.screen import Screen, ModalScreen
from textual import events
from textual.reactive import reactive, var
from textual.binding import Binding
from textual.message import Message

from rich.json import JSON
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

# Import the comparison functionality
from post import PipelineMetricsCollection, SampleMetrics
import re

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
        ("j", "cursor_down", "Down"),
        ("k", "cursor_up", "Up"),
        ("enter", "select_cursor", "Select"),
        ("space", "select_cursor", "Select"),
    ]
    
    def __init__(self, path: Path, **kwargs):
        super().__init__(path, **kwargs)
        self.selected_samples = set()
        self.structure_type = self._detect_structure()
        self.sample_nodes = {}  # Track sample nodes for styling
    
    def _detect_structure(self) -> str:
        """Detect if this is 'workdir' or 'experiment' structure."""
        try:
            # Check if immediate subdirectories contain pipeline_summary.json
            direct_samples = any(
                (subdir / "pipeline_summary.json").exists() 
                for subdir in self.path.iterdir() 
                if subdir.is_dir() and not subdir.name.startswith('.')
            )
            
            if direct_samples:
                return "experiment"  # Direct samples: experiment/sample/pipeline_summary.json
            else:
                return "workdir"     # Nested: workdir/experiment/sample/pipeline_summary.json
        except PermissionError:
            return "workdir"  # Default fallback
    
    def filter_paths(self, paths: Iterable[Path]) -> Iterable[Path]:
        """Filter based on detected structure."""
        for path in paths:
            if not path.is_dir() or path.name.startswith('.'):
                continue
                
            if self.structure_type == "experiment":
                # Show sample directories with pipeline_summary.json
                if (path / "pipeline_summary.json").exists():
                    yield path
            else:  # workdir structure
                # Show experiment directories that contain valid samples
                if self._is_experiment_dir(path):
                    yield path
                elif self._is_valid_sample_dir(path):
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
        if self.structure_type == "experiment":
            # Just the sample name
            return sample_path.name
        else:
            # experiment/sample format
            return f"{sample_path.parent.name}/{sample_path.name}"
    
    def get_selected_sample_paths(self) -> list[Path]:
        """Get paths to all selected samples."""
        selected_paths = []
        for sample_id in self.selected_samples:
            if self.structure_type == "experiment":
                sample_path = self.path / sample_id
            else:
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
        if self.structure_type == "experiment":
            sample_path = self.path / sample_id
        else:
            if '/' not in sample_id:
                return False
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
        ("q", "return_to_main", "Return to Main"),
        ("r", "run_comparison", "Run Comparison"),
    ]

    def __init__(self, collection: PipelineMetricsCollection, preselected_samples=None):
        super().__init__()
        self.collection = collection
        self.selected_samples = set(preselected_samples or [])
        self.selected_samples = set(self._normalize_ids())

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

                # Create SelectionList with preselected samples
                selection_options = [(sample_id, sample_id) for sample_id in self.selected_samples]
                selection_list = SelectionList[str](*selection_options, id="sample-selection-list")
                
                # Pre-select all items
                for sample_id in self.selected_samples:
                    selection_list.select(sample_id)
                
                yield selection_list

            with Vertical(id="metric-buttons"):
                yield Button("Compare Valid UMI %", classes="metric-button", id="compare-valid-umi")
                yield Button("Compare Read Coverage", classes="metric-button", id="compare-read-coverage")

            # Results display
        with Vertical(id="right-panel"):
            yield ComparisonDataTable(id="comparison-table")

        yield Footer()

    def _run_comparison(self) -> None:
        """Run the selected comparison."""
        # Get currently selected samples from SelectionList
        selection_list = self.query_one("#sample-selection-list", SelectionList)
        selected_indices = selection_list.selected
        
        if not selected_indices:
            self.notify("Please select at least one sample to compare", severity="warning")
            return

        # Get the sample IDs for selected indices
        all_samples = list(self.selected_samples)
        current_selected = {all_samples[i] for i in selected_indices}

        # Default to valid UMI comparison
        df = self.collection.get_valid_umi_stats()
        df = df.filter(pl.col("sample_id").is_in(list(current_selected)))

        if df.height == 0:
            self.notify("No valid data for selected samples", severity="warning")
            return

        self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Valid UMI Comparison")

    def _normalize_ids(self) -> List:
        """Normalizes ids to handle paths at different levels"""
        normalized_ids = []
        for sample_id in self.selected_samples:
            if '/' in sample_id:
                normalized_ids.append(sample_id.split('/')[-1])
            else:
                normalized_ids.append(sample_id)
        return normalized_ids

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
        ("q", "quit", "Quit"),
        ("r", "refresh", "Refresh"),
        ("f", "view_figures", "View Figures"),
        ("c", "compare_samples", "Compare Samples"),
        ("e", "toggle_files", "Toggle Files"),
        ("a", "select_all", "Select All Files"),
        ("d", "deselect_all", "Clear Selection"),
    ]
    show_tree = var(False)


    def action_deselect_all(self) -> None:
        """Select All Files"""
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        tree.deselect_all_samples()

    def action_select_all(self) -> None:
        """Select All Files"""
        tree = self.query_one("#experiment-tree", SmartExperimentDirectoryTree)
        selected_n = tree.select_all_samples()
        self.notify(f"Selected {selected_n} samples")

    def action_toggle_files(self) -> None:
        """Called in response to key binding."""
        self.show_tree = not self.show_tree
    
    def watch_show_tree(self, show_tree: bool) -> None:
        """Called when show_tree is modified."""
        self.set_class(show_tree, "-show-raw")

    def __init__(self, experiment_dir: str):
        """Initialize the app with the experiment directory."""
        super().__init__()
        self.experiment_dir = Path(experiment_dir).resolve()
        if not self.experiment_dir.exists() or not self.experiment_dir.is_dir():
            raise ValueError(f"Experiment directory {experiment_dir} does not exist.")
        self.current_figures = []
        self.current_sample_dir = None

        # For sample comparison
        self.selected_samples = set()
        self.metrics_collection = PipelineMetricsCollection()
        self.pipeline_files = []

    def compose(self) -> ComposeResult:
        """Compose the app layout."""
        
        #yield Static(title_art, id="title-art")
        yield Header(show_clock=True)
        
        with Vertical(id='experiments'):
            yield SmartExperimentDirectoryTree(self.experiment_dir, id="experiment-tree")

        with Vertical(id="metrics-container"):
            yield MetricsTable(id="all-metrics", stage="all")
        
            # Add a JSONDisplay to show the full JSON data
            yield JSONDisplay(id="raw-json-display")

            with Horizontal(id="control-buttons"):
                yield Button("Select All", id="select-all-button")
                yield Button("Clear Selection", id="deselect-all-button")
                yield Button("Compare Samples", id="compare-samples-button")
                yield Button("Select Sample", id="select-sample-button")
                yield Button("View Figures", id="figures-button")

        yield Footer()

    def on_mount(self) -> None:
        """Setup when the app is mounted."""
        # Remove the SelectionList references since we're not using it anymore
        self.theme = "gruvbox"

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
    """Run the app."""
    if len(sys.argv) > 1:
        experiment_dir = sys.argv[1]
    else:
        experiment_dir = os.environ.get("FRACTURE_EXPERIMENT_DIR", 
                                        os.path.expanduser("~/src/fracture-app/local-toy"))
    
    app = FractureExplorer(experiment_dir)
    app.run()


if __name__ == "__main__":
    main()
