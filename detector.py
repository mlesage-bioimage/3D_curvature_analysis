"""
DETECTOR - Same logic as your original detect_stage_from_filename
"""

from config import GROUPS_ORDERED, GROUP_PATTERNS


def detect_group_from_filename(filename):
    """
    Detect group from filename - SAME LOGIC as your original detect_stage_from_filename
    """
    filename_lower = filename.lower()
    
    # Check each group's patterns
    for group in GROUPS_ORDERED:
        for pattern in GROUP_PATTERNS[group]:
            if pattern.lower() in filename_lower:
                return group
    
    return "unknown"


def get_groups_ordered():
    """Return groups in order"""
    return GROUPS_ORDERED


# Class for compatibility with your existing code structure
class StageDetector:
    def __init__(self):
        self.order = GROUPS_ORDERED
    
    def detect(self, filename):
        return detect_group_from_filename(filename)
    
    def get_order(self):
        return self.order