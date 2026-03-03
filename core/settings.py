import json
import os

SETTINGS_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "settings.json")

DEFAULTS = {
    "min_overlap": 3,
    "max_overlap": 10,
    "max_mismatches": 1,
    "consider_ambiguity": False,
    "high_risk_min_overlap": 4,
    "high_risk_max_mismatches": 0,
    "medium_risk_min_overlap": 2,
    "medium_risk_max_mismatches": 1,
}


def load_settings():
    """Load settings from settings.json. Creates file with defaults if missing."""
    if os.path.exists(SETTINGS_FILE):
        with open(SETTINGS_FILE, "r") as f:
            saved = json.load(f)
        # Merge with defaults so new keys are always present
        merged = {**DEFAULTS, **saved}
        return merged
    else:
        save_settings(DEFAULTS)
        return dict(DEFAULTS)


def save_settings(data):
    """Save settings dict to settings.json."""
    with open(SETTINGS_FILE, "w") as f:
        json.dump(data, f, indent=4)
