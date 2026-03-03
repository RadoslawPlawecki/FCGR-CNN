"""
@author: Radosław Pławecki
"""

import numpy as np
import re
import matplotlib.pyplot as plt
import math

class FCGR:
    """
    Class to implement the Frequency Chaos Game Representation (FCGR) to encode DNA sequences. 
    """
    start = (0.5, 0.5)
    ade = (0, 1)
    cyt = (1, 1)
    thy = (0, 0)
    gua = (1, 0)

    def __init__(self, sequence, k_mer):
        """
        Method to initialize a class.
        
        :param self: a reference.
        :param sequence: DNA sequence.
        :param k_mer: length of DNA pieces.
        """
        if self.__is_valid_dna(sequence=sequence) is False:
            raise ValueError("Sequence is not correct! Allowed characters: A, C, T, G.")
        self.sequence = sequence.upper()
        if k_mer < 2 or k_mer > 6:
            raise ValueError("k-mer takes values from the range [3,6]!")
        if k_mer > len(sequence):
            raise ValueError("k-mer cannot be longer than a sequence!")
        self.k_mer = k_mer
        self.img_size = int(np.sqrt(4 ** self.k_mer))

    @staticmethod
    def __is_valid_dna(sequence):
        """
        Method to check if a DNA sequence is valid.
        
        :param sequence: DNA sequence.
        :return bool: True, if a DNA sequence is valid, otherwise - False.
        """
        return bool(re.fullmatch(r"[ACGTacgt]+", sequence))

    @staticmethod
    def __movement_rule(actual, fixed):
        """
        Method to calculate new point.
        
        :param actual: actual point.
        :param fixed: point to move.
        :return new point:
        """
        return ((actual[0] + fixed[0]) / 2, (actual[1] + fixed[1]) / 2)

    def compute_points(self):
        """
        Method to compute all the points.
        
        :param self: a reference.
        :return points: computed points.
        """
        seq = self.sequence
        point = self.start
        points = []
        for char in seq:
            if char == 'A':
                new_point = self.__movement_rule(point, self.ade)
            elif char == 'C':
                new_point = self.__movement_rule(point, self.cyt)
            elif char == 'G':
                new_point = self.__movement_rule(point, self.gua)
            elif char == 'T':
                new_point = self.__movement_rule(point, self.thy)
            point = new_point
            points.append(point)
        return points

    def __init_matrix(self):
        """
        Method to initialize a matrix.
        
        :param self: a reference.
        :return initialized matrix with zeros:
        """
        rows = cols = self.img_size
        return np.zeros((rows, cols))
    
    def point_pixel_map(self):
        """
        Method to map computed points to pixels.
        
        :param self: a reference.
        :return locations of pixels:
        """
        points = self.compute_points()
        const = self.img_size
        pixels = []
        for i in range(len(points)):
            if i < self.k_mer - 1:
                continue
            pixel = (math.floor(points[i][1] * const), math.floor(points[i][0] * const))
            pixels.append(pixel)
        return pixels

    def fill_matrix(self):
        """
        Method to fill a matrix with values.
        
        :param self: a reference.
        :return filled matrix:
        """
        matrix = self.__init_matrix()
        pixels = self.point_pixel_map()
        for (r, c) in pixels:
            matrix[r][c] += 1
        return matrix
    
    def display_matrix(self):
        """
        Method to display a matrix.
        
        :param self: a reference.
        """
        matrix = np.array(self.fill_matrix(), dtype=np.float32)
        plt.figure()
        plt.imshow(matrix, cmap='gray')
        plt.colorbar()
        plt.tight_layout()
        plt.show()
