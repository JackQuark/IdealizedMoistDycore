from dycore import Any

import tomllib

class Config(object):
    """Handler for dycore configuration files in TOML format."""
    @classmethod
    def from_toml(cls, filepath: str):
        with open(filepath, "rb") as f:
            config_dict = tomllib.load(f)
        return cls(config_dict)
    
    def __init__(self, config_dict: dict[str, dict[str, Any]] = None):
        super(Config, self).__init__()
        if config_dict is None:
            config_dict = {}
        self.config = config_dict
        
    def output_toml(self, filepath: str):
        with open(filepath, "wb") as f:
            toml_bytes = tomllib.dumps(self.config).encode("utf-8")
            f.write(toml_bytes)
