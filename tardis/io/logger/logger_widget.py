import logging
import re
import pandas as pd
from tardis.util.environment import Environment
import panel as pn
pn.extension(comms="ipywidgets")


def create_logger_columns(start_height=10, max_height=300):
    """Create logger scroll columns with dynamic height.
    
    Returns
    -------
    dict
        Dictionary of scroll columns for each log level.
    """
    # Create scroll columns for each log level
    columns = {}
    
    for level in ["INFO", "WARNING/ERROR", "DEBUG"]:
        column = pn.pane.HTML(
            "",
            height=10,  # Start small
            styles={
                'border': '1px solid #ddd',
                'background-color': 'white',
                'overflow-y': 'auto',
                'overflow-x': 'auto',
                'font-family': 'monospace',
                'padding': '8px',
                'white-space': 'pre-wrap'
            },
            sizing_mode='stretch_width'
        )
        column.max_log_entries = 1000
        column._start_height = start_height
        column._max_height = max_height
        column._log_content = ""  # Store content as string
        columns[level] = column
    
    return columns


class PanelWidgetLogHandler(logging.Handler):
    """Log handler for logging to scroll columns.
    
    Parameters
    ----------
    log_columns : dict
        Dictionary of scroll columns for each log level.
    colors : dict
        Dictionary mapping log levels to display colors.
    display_widget : bool, optional
        Whether to display logs in the widget. Defaults to True.
    display_handles : dict, optional
        Dictionary of display handles for each column (jupyter environment).
    """
    def __init__(self, log_columns, colors, display_widget=True, display_handles=None):
        super().__init__()
        self.log_columns = log_columns
        self.colors = colors
        self.display_widget = display_widget
        self.display_handles = display_handles or {}
        self.environment = Environment.get_current_environment()
        
        self.stream_handler = None
        if not self.display_widget:
            self.stream_handler = logging.StreamHandler()
            self.stream_handler.setFormatter(logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"))

    def emit(self, record):
        """Process and emit a log record.
        
        Parameters
        ----------
        record : logging.LogRecord
            The log record to process and display.
        """
        if not self.display_widget:
            if self.stream_handler:
                self.stream_handler.emit(record)
            return

        if isinstance(record.msg, pd.DataFrame):
            df_str = record.msg.to_string(index=True)
            html_output = f'<pre style="white-space: pre; margin: 0;">{df_str}</pre>'
        else:
            log_entry = self.format(record)
            clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
            html_output = self._format_html_output(clean_log_entry, record)

        self._emit_to_columns(record.levelno, html_output)

    @staticmethod
    def _remove_ansi_escape_sequences(text):
        """Remove ANSI escape sequences from string.
        
        Parameters
        ----------
        text : str
            The text containing ANSI escape sequences.
            
        Returns
        -------
        str
            Cleaned text with ANSI escape sequences removed.
        """
        ansi_escape = re.compile(r"\x1B[@-_][0-?]*[ -/]*[@-~]")
        return ansi_escape.sub("", text)

    def _format_html_output(self, log_entry, record):
        """Format log entry as HTML with appropriate styling.
        
        Parameters
        ----------
        log_entry : str
            The log entry text to format.
        record : logging.LogRecord
            The log record containing level information.
            
        Returns
        -------
        str
            HTML-formatted log entry.
        """
        color = self.colors.get(record.levelno, self.colors["default"])
        parts = log_entry.split(" ", 2)
        if len(parts) > 2:
            prefix, levelname, message = parts
            return f'<span>{prefix}</span> <span style="color: {color}; font-weight: bold;">{levelname}</span> {message}'
        return log_entry

    def _emit_to_columns(self, level, html_output):
        """Add log entry to appropriate scroll columns using string concatenation.
        
        Parameters
        ----------
        level : int
            The logging level.
        html_output : str
            The HTML-formatted log message.
        """
        level_to_output = {
            logging.WARNING: "WARNING/ERROR",
            logging.ERROR: "WARNING/ERROR",
            logging.CRITICAL: "WARNING/ERROR",
            logging.INFO: "INFO",
            logging.DEBUG: "DEBUG"
        }

        # Add to specific level column using string concatenation
        output_key = level_to_output.get(level)
        if output_key and output_key in self.log_columns:
            level_column = self.log_columns[output_key]
            
            # Wrap in div and concatenate to existing content
            html_wrapped = f"<div style='margin: 2px 0; padding: 2px 0;'>{html_output}</div>"
            current_content = getattr(level_column, '_log_content', '') or level_column.object or ''
            new_content = current_content + html_wrapped if current_content else html_wrapped
            
            # Store content and update column
            level_column._log_content = new_content
            level_column.object = new_content
            
            # Dynamic height adjustment
            self._adjust_column_height(level_column, new_content)
            
            # Trim old entries if needed (count divs)
            if new_content.count('<div') > level_column.max_log_entries:
                # Keep only the last max_log_entries
                divs = new_content.split('<div')
                trimmed = '<div'.join(divs[-level_column.max_log_entries:])
                level_column._log_content = trimmed
                level_column.object = trimmed
            
            # Update display handle in jupyter environment
            if (self.environment == 'jupyter' and output_key in self.display_handles 
                and self.display_handles[output_key] is not None):
                self.display_handles[output_key].update(level_column)
    
    def _adjust_column_height(self, column, content=None):
        """Dynamically adjust column height based on content.
        
        Parameters
        ----------
        column : panel.pane.HTML
            The column to adjust.
        content : str, optional
            The content to count entries from. If None, uses column._log_content.
        """
        if hasattr(column, '_start_height') and hasattr(column, '_max_height'):
            # Count entries by counting divs in content
            if content is None:
                content = getattr(column, '_log_content', '') or column.object or ''
            entry_count = content.count('<div') if content else 0
            
            # Estimate height needed: ~25px per log entry
            estimated_height = max(column._start_height, entry_count * 25)
            # Cap at max height
            new_height = min(estimated_height, column._max_height)
            column.height = new_height
    
    def close(self):
        """Close the log handler.
        """
        super().close()
        if self.stream_handler:
            self.stream_handler.flush()
            self.stream_handler.close()
            self.stream_handler = None


