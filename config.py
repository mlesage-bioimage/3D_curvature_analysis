"""
CONFIGURATION FILE - Change only this file for different experiments
"""

# ========== EXPERIMENT SETUP ==========
EXPERIMENT_NAME = "Nucleus_Morphology_Analysis"

# ========== INPUT/OUTPUT SETTINGS ==========
INPUT_FOLDER = r"path/to/your/images"
OUTPUT_EXCEL = r"path/to/your/results.xlsx"

# ========== GROUP DETECTION SETTINGS ==========
# List your groups in order (for Excel sheets)
GROUPS_ORDERED = ["DMSO", "DRB"]

# Patterns to detect each group (order matters - longer first if needed)
# Each pattern will be looked for in the filename
GROUP_PATTERNS = {
    "DMSO": ["DMSO"],      # Files containing "DMSO" go to DMSO group
    "DRB": ["DRB"],        # Files containing "DRB" go to DRB group
}

# ========== MICROSCOPE SETTINGS ==========
ORIGINAL_SPACING = (0.258, 0.258, 1.0)  # Your original spacing
TARGET_SPACING = 0.258
CURVATURE_RADIUS = 1.0
CURVATURE_BINS = 100
Z_EXCLUSION = 0.07

# ========== PROCESSING SETTINGS ==========
USE_PARALLEL = True
MAX_WORKERS = 2



# =============================================================================
# For embryo stages :
# GROUPS_ORDERED = ["zygote", "2cells", "4cells", "8cells", "16cells", "32cells", "blastocyst"]
# GROUP_PATTERNS = {
#     "zygote": ["zygote"],
#     "2cells": ["2cells"],
#     "4cells": ["4cells"],
#     "8cells": ["8cells"],
#     "16cells": ["16cells"],
#     "32cells": ["32cells"],
#     "blastocyst": ["blasto"],
# }
# =============================================================================

# =============================================================================
# For different treatments :
# GROUPS_ORDERED = ["Control", "DrugA", "DrugB"]
# GROUP_PATTERNS = {
#     "Control": ["control", "DMSO", "vehicle"],
#     "DrugA": ["drugA", "A_treated"],
#     "DrugB": ["drugB", "B_treated"],
# }
# =============================================================================
