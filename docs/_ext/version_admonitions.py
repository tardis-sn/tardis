"""
Custom version mark admonition directives for TARDIS documentation.
Provides consistent admonition boxes for RTX and Classic architectures.
"""
from docutils.parsers.rst import Directive
from docutils import nodes
from pathlib import Path
import toml


# Load version mark texts from external file
_config_path = Path(__file__).parent / "version_marks.toml"
_VERSION_MARKS = toml.load(_config_path)


class TardisVersionMark(Directive):
    """
    Base class for TARDIS version-specific admonition directives.
    
    Subclasses should define the `version` class attribute.
    """
    
    has_content = False
    required_arguments = 0
    optional_arguments = 0
    version = None  # Override in subclasses
    
    def run(self):
        """Generate the admonition node."""
        if self.version not in _VERSION_MARKS:
            raise ValueError(f"Unknown version mark: {self.version}")
        
        config = _VERSION_MARKS[self.version]
        
        # Create admonition node with version-specific class
        admonition_node = nodes.admonition()
        admonition_node['classes'].append(self.version)
        admonition_node['classes'].append('admonition')
        
        # Add title
        title_text = config['title']
        title = nodes.title(title_text, title_text)
        admonition_node += title
        
        # Add content paragraph
        para = nodes.paragraph()
        para += nodes.Text(config['text'])
        admonition_node += para
        
        return [admonition_node]


class RTXAdmonition(TardisVersionMark):
    """
    RTX architecture admonition directive.
    
    Usage in RST files:
        .. rtx-admonition::
    """
    version = "rtx"


class ClassicAdmonition(TardisVersionMark):
    """
    Classic architecture admonition directive.
    
    Usage in RST files:
        .. classic-admonition::
    """
    version = "classic"


def setup(app):
    """Setup the version mark admonition directives."""
    app.add_directive("rtx-admonition", RTXAdmonition)
    app.add_directive("classic-admonition", ClassicAdmonition)
    
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
