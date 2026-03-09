"""
@author: Radosław Pławecki
"""

from chunks_preprocessor import ChunksPreprocessor

input_path = "data/clinvar/chunks"
output_path = "data/clinvar/03_preprocessed_clinvar.csv"

cp = ChunksPreprocessor(dir_path=input_path)
cp.preprocess_chunks(output_path=output_path)
