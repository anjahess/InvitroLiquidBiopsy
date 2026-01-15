"""

Constants for In vitro Liquid Biopsy paper

Author: Anja Hess
Date: 2026 JAN 15


"""
import os

script_dir = os.path.dirname(__file__)
PACKAGE_DIR = script_dir.rsplit('scripts/utils')[0]
SOURCE_DATA_DIR = f"{PACKAGE_DIR}/sourcedata/"
UTILS_DIR = f"{PACKAGE_DIR}/scripts/utils/"

RESULTS_DIR = f"{PACKAGE_DIR}results/"
FIGURE_DIR = f"{RESULTS_DIR}figures/"
TABLE_DIR = f"{RESULTS_DIR}tables/"
