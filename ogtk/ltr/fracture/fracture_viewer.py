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
from rich.json import JSON
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

# Import the comparison functionality
from fracture_compare import PipelineMetricsCollection, SampleMetrics

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

    def update_metrics(self, metrics: Dict[str, Any], stage: str) -> None:
        """Update the metrics data."""
        self.metrics = metrics
        self.stage = stage
        self.update(self._render_metrics())

    def _format_number(self, value, metric_name) -> str:
        """Format numeric values with underscore separators for thousands."""
        if isinstance(value, (int, float)):
            # For percentages, rates, and similar metrics, format with 2 decimal places
            if isinstance(value, float) and any([i in metric_name for i in ["fraction", "rate", "mean", "median"]]):
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

    def _render_metrics(self) -> Panel:
        """Render the metrics as a table."""
        if not self.metrics:
            return Panel("No metrics selected")

        table = Table(title=f"{self.stage.capitalize()} Metrics")
        table.add_column("Metric")
        table.add_column("Value")

        for key, value in self.metrics.items():
            formatted_value = self._format_number(value, key)
            table.add_row(key, formatted_value)

        return Panel(table, title=f"{self.stage.capitalize()} Metrics")


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
            table.add_row(*[str(row[col]) for col in self.df.columns])

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

    def action_return_to_main(self) -> None:
        """Return to main screen."""
        self.app.pop_screen()

    def action_run_comparison(self) -> None:
        """Run the comparison with selected samples."""
        self._run_comparison()

    def compose(self) -> ComposeResult:
        """Compose the screen layout."""
        with Container(id="comparison-container"):
            yield Header(show_clock=True)

            with Container(id="sample-selection"):
                yield Label("Pre-selected samples to compare:", id="selection-label")

                pre_selected = [(sample, sample, True) for sample in self.selected_samples]
                yield SelectionList(*pre_selected, id="pre-selected")

                yield Button("Compare Selected Samples", id="compare-button")

            with Container(id="comparison-metrics"):
                # Add preset comparison buttons
                with Horizontal():
                    yield Button("Compare Valid UMI %", classes="metric-button", id="compare-valid-umi")
                    yield Button("Compare Read Coverage", classes="metric-button", id="compare-read-coverage")

                # Results display
                yield ComparisonDataTable(id="comparison-table")

            yield Footer()

    def on_mount(self) -> None:
        """Setup the screen when mounted."""
        self._update_sample_list()

    def _update_sample_list(self) -> None:
        """Update the sample selection grid."""
        samples = self.selected_samples


        # Calculate the grid layout - we want roughly a square grid
        import math
        cols = max(1, int(math.sqrt(len(samples))))


        # Add checkboxes for each sample
        for sample_id in samples:
            sample_label = f"{sample_id}"
            checkbox = Checkbox(sample_label, classes="sample-checkbox", id=f"sample-{sample_id}")
            checkbox.value = True
            checkbox.add_class("checkbox-selected")


    def _run_comparison(self) -> None:
        """Run the selected comparison."""
        if not self.selected_samples:
            self.notify("Please select at least one sample to compare", severity="warning")
            return

        # Default to valid UMI comparison
        df = self.collection.get_valid_umi_stats()

        # Filter to selected samples only
        df = df.filter(pl.col("sample_id").is_in(list(self.selected_samples)))

        if df.height == 0:
            self.notify("No valid data for selected samples", severity="warning")
            return

        self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Valid UMI Comparison")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        button_id = event.button.id

        if button_id == "compare-button":
            self._run_comparison()

        elif button_id == "compare-valid-umi":
            if not self.selected_samples:
                self.notify("Please select at least one sample to compare", severity="warning")
                return

            df = self.collection.get_valid_umi_stats()
            df = df.filter(pl.col("sample_id").is_in(list(self.selected_samples)))

            if df.height == 0:
                self.notify("No valid UMI data for selected samples", severity="warning")
                return

            self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Valid UMI Comparison")

        elif button_id == "compare-read-coverage":
            if not self.selected_samples:
                self.notify("Please select at least one sample to compare", severity="warning")
                return

            df = self.collection.calculate_read_coverage()
            df = df.filter(pl.col("sample_id").is_in(list(self.selected_samples)))

            if df.height == 0:
                self.notify("No read coverage data for selected samples", severity="warning")
                return

            self.query_one("#comparison-table", ComparisonDataTable).update_data(df, "Read Coverage Comparison")

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        """Handle checkbox selection changes."""
        checkbox_id = event.checkbox.id
        if checkbox_id and checkbox_id.startswith("sample-"):
            sample_id = checkbox_id[7:]  # Remove "sample-" prefix

            if event.value:
                self.selected_samples.add(sample_id)
            else:
                self.selected_samples.discard(sample_id)


class FractureExplorer(App):
    """A Textual app to explore Fracture pipeline outputs."""

    CSS_PATH = "fracture_viewer.tcss"

    BINDINGS = [
        ("q", "quit", "Quit"),
        ("r", "refresh", "Refresh"),
        ("f", "view_figures", "View Figures"),
        ("c", "compare_samples", "Compare Samples"),
        ("e", "toggle_files", "Toggle Files"),
    ]
    show_tree = var(False)

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
        yield Header(show_clock=True)
        #yield FilteredDirectoryTree(self.experiment_dir, id="dir-tree")
        yield JSONDisplay(id="raw-json-display")

        with Vertical(id='experiments'):
            yield ExperimentList(id="experiment-list", path=self.experiment_dir)

        with Vertical(id="metrics-container"):
            yield MetricsTable(id="all-metrics", stage="all")

            # Add a JSONDisplay to show the full JSON data
            #yield JSONDisplay(id="raw-json-display")

            with Horizontal(id="control-buttons"):
                yield Button("Compare Samples", id="compare-samples-button")
                yield Button("Select Sample", id="select-sample-button")
                yield Button("View Figures", id="figures-button")

        yield Footer()

    def on_mount(self) -> None:
        self.query_one(SelectionList).border_title = "Select samples to compare"
        #self.query_one(Pretty).border_title = "Selected games"

    @on(events.Mount)
    @on(SelectionList.SelectedChanged)
    def update_selected_view(self) -> None:
        selected_set = self.query_one(SelectionList).selected 
        #self.query_one(Pretty).update(selected_set)
        self.selected_samples = selected_set

    @on(SelectionList.SelectionHighlighted)
    def on_experiment_highlighted(self, event:SelectionList.SelectionHighlighted):
        experiment_list = self.query_one(SelectionList)
        highlighted_index = event.selection.value
        #self.query_one(Pretty).update(highlighted_index)

    def on_selection_list_selection_highlighted(self, event: SelectionList.SelectionHighlighted) -> None:
        """Load and display JSON data when navigating with keyboard in the selection list."""
        sample_name = event.selection.value
        
        # Find the sample directory
        sample_dir = None
        for path in self.experiment_dir.iterdir():
            if path.is_dir() and path.name == sample_name and (path/'pipeline_summary.json').exists():
                sample_dir = path
                break
        
        if sample_dir:
            # Load the JSON data
            try:
                with open(sample_dir / "pipeline_summary.json", "r") as f:
                    data = json.load(f)
                
                # Update the current sample directory
                self.current_sample_dir = sample_dir
                
                # Reset figures list
                self.current_figures = []
                figures_button = self.query_one("#figures-button", Button)
                figures_button.remove_class("has-figures")
                
                # Update selection button
                sample_id = sample_dir.name
                select_button = self.query_one("#select-sample-button", Button)

                if sample_id in self.selected_samples:
                    select_button.label = f"Unselect '{sample_id}'"
                else:
                    select_button.label = f"Select '{sample_id}'"
                select_button.display = True
                
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
                        all_metrics, "all"
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
        tree = self.query_one(FilteredDirectoryTree)
        tree.reload()

    def action_view_figures(self) -> None:
        """View figures for the current sample."""
        if self.current_figures and self.current_sample_dir:
            self.push_screen(FigureViewer(self.current_figures, self.current_sample_dir))
        else:
            self.notify("No figures available for the current sample", severity="warning")

    def action_compare_samples(self) -> None:
        """Open the sample comparison screen."""
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

            self.push_screen(ComparisonScreen(self.metrics_collection, preselected_samples=self.selected_samples))
        except Exception as e:
            self.notify(f"Error loading metrics for comparison: {str(e)}", severity="error")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button press events."""
        button_id = event.button.id

        if button_id == "figures-button":
            self.action_view_figures()

        elif button_id == "compare-samples-button":
            self.action_compare_samples()

        elif button_id == "select-sample-button":
            if self.current_sample_dir:
                sample_id = self.current_sample_dir.name
                select_button = self.query_one("#select-sample-button", Button)

                if sample_id in self.selected_samples:
                    # Unselect the sample
                    self.selected_samples.remove(sample_id)
                    self.notify(f"Sample '{sample_id}' removed from selection")
                    select_button.label = f"Select '{sample_id}'"

                    # Update directory tree styling
                    try:
                        # Find and update the tree node styling
                        tree = self.query_one(FilteredDirectoryTree)
                        for node in tree.nodes.values():
                            if node.data.path == self.current_sample_dir:
                                node.remove_class("selected-sample")
                                break
                    except Exception:
                        pass  # Ignore if we can't find the node
                else:
                    # Select the sample
                    self.selected_samples.add(sample_id)
                    self.notify(f"Sample '{sample_id}' added to selection")
                    select_button.label = f"Unselect '{sample_id}'"

                    # Update directory tree styling
                    try:
                        # Find and update the tree node styling
                        tree = self.query_one(FilteredDirectoryTree)
                        for node in tree.nodes.values():
                            if node.data.path == self.current_sample_dir:
                                node.add_class("selected-sample")
                                break
                    except Exception:
                        pass  # Ignore if we can't find the node
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
