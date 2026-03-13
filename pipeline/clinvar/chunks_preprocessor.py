"""
@author: Radosław Pławecki
"""

import os
import re
import pandas as pd
from tqdm import tqdm


class ChunksPreprocessor:
    """
    Class to process the chunks of ClinVar data.
    """
    _COLS = ["GRCh38Chromosome", "GRCh38Location", "label", "GenomicReference"]
    _CHR_MAP = {"X": 23, "MT": 0}
    _HGVS_PATTERN = re.compile(
        r'^(?P<accession>[^:]+):'
        r'(?P<coord_type>[gm])\.'
        r'(?P<position>\d+)'
        r'(?P<ref>[ACGT])>(?P<mut>[ACGT])$')

    def __init__(self, dir_path):
        """
        Initialize the class instance.
        """
        if not os.path.isdir(dir_path):
            raise ValueError("The path should be a directory with data chunks!")
        self.dir_path = dir_path

    def _load(self, file_path) -> pd.DataFrame:
        """
        Load selected columns from a CSV file into a DataFrame.
        """
        return pd.read_csv(file_path, delimiter=';', usecols=self._COLS, dtype={"GenomicReference": "string"})

    def _parse_gen_ref(self, variant: str):
        """
        Parse genomic references and split them into their components.
        """
        match = self._HGVS_PATTERN.match(variant)
        if not match:
            raise ValueError(f"Invalid variant: '{variant}'!")
        return match.groupdict()
    
    def _clean(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean data based on genomic references, e.g. missing values.
        """
        return df[df["GenomicReference"].str.match(self._HGVS_PATTERN)]
    
    def _remove_mismatches(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove mismatches between GRCh38Location and position in genomic references.
        """
        pos = df['GenomicReference'].str.extract(r'(?:g|m)\.(\d+)')[0]
        pos = pd.to_numeric(pos, errors='coerce').astype('Int64')
        grch38_loc = pd.to_numeric(df["GRCh38Location"], errors='coerce').astype('Int64')
        return df.loc[grch38_loc == pos]

    def _expand_gen_ref(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add new columns to the DataFrame based on genomic reference components.
        """
        df_parsed = df["GenomicReference"].apply(self._parse_gen_ref).apply(pd.Series)
        return pd.concat([df.reset_index(drop=True), df_parsed.reset_index(drop=True)], axis=1)

    def _map_chromosomes(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Map chromosome string values to numeric codes (X → 23, MT → 0).
        """
        df = df.copy()
        df["GRCh38Chromosome"] = df["GRCh38Chromosome"].replace(self._CHR_MAP)
        return df

    def _organise_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Rename and reorganise columns of the DataFrame.
        """
        df.rename(columns={"GRCh38Chromosome": "chr", "GRCh38Location": "loc"}, inplace=True)
        df = df[["accession", "chr", "loc", "ref", "mut", "label"]]
        return df

    def _process_chunk(self, file_path) -> pd.DataFrame:
        """
        Process chunks using pipeline.
        """
        df = self._load(file_path)
        df = self._clean(df)
        df = self._remove_mismatches(df)
        if df.empty:
            return df
        df = self._expand_gen_ref(df)
        if df.empty:
            return df
        df = self._map_chromosomes(df)
        df = self._organise_data(df)
        return df

    def preprocess_chunks(self, output_path=None):
        """
        Execute the pipeline and save the output to a CSV file.
        """
        first = True
        print("[INFO] Processing chunks...")
        for chunk_path in tqdm(os.listdir(self.dir_path)):
            file_path = os.path.join(self.dir_path, chunk_path)
            if not os.path.isfile(file_path):
                continue
            df_chunk = self._process_chunk(file_path)
            if output_path:
                df_chunk.to_csv(output_path, mode='w' if first else 'a', header=first, index=False, sep=';')
            first = False
        print(f"[INFO] Data was exported!")
