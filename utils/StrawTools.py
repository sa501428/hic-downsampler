import numpy as np
import strawC
from scipy import sparse
from scipy.sparse import coo_matrix


def extract_all_data_hic(url, normalization, resolution, chroms1, chroms2, chrom_dot_sizes):
    overall_array = None
    for chr1 in chroms1:
        current_row = None
        for chr2 in chroms2:
            need_to_flip = False
            achr1, achr2 = chr1, chr2

            if int(chr2) < int(chr1):
                need_to_flip = True
                achr1, achr2 = chr2, chr1

            max_rows = int(chrom_dot_sizes[achr1] / resolution) + 1
            max_cols = int(chrom_dot_sizes[achr2] / resolution) + 1

            result = strawC.strawC("observed", normalization, url, achr1, achr2, 'BP', resolution)

            result_x = list()
            result_y = list()
            result_val = list()

            for i in range(len(result)):
                result_x.append(int(result[i].binX))
                result_y.append(int(result[i].binY))
                result_val.append(result[i].counts)
            if int(chr2) == int(chr1):
                for i in range(len(result)):
                    if int(result[i].binX) != int(result[i].binY):
                        result_x.append(int(result[i].binY))
                        result_y.append(int(result[i].binX))
                        result_val.append(result[i].counts)

            row_indices = np.asarray(result_x) / resolution
            col_indices = np.asarray(result_y) / resolution
            data = result_val
            mat_coo = coo_matrix((data, (row_indices.astype(int), col_indices.astype(int))),
                                 shape=(max_rows, max_cols))

            # put data into a dense matrix format
            dense = mat_coo.toarray()
            dense[np.isnan(dense)] = 0  # set NaNs to zero

            if need_to_flip:
                dense = dense.T
            print('.')

            if current_row is None:
                current_row = dense
            else:
                current_row = np.hstack((current_row, dense))

        if overall_array is None:
            overall_array = current_row
        else:
            overall_array = np.vstack((overall_array, current_row))
    return overall_array


def get_chromosome_dot_sizes(hic_file: str):
    cds = {}
    chromosome_dot_sizes = strawC.getChromosomes(hic_file)
    for chromosome in chromosome_dot_sizes:
        temp_name = str(chromosome.name).lower()
        if "all" in temp_name or "mt" in temp_name or "y" in temp_name:
            continue
        try:
            chromosome_index = int(temp_name.replace("chr", ""))
            cds[chromosome_index] = chromosome.length
        except:
            pass
    return cds


class HiCFile:
    def __init__(self, filepath: str, resolutions: list, norm: str):
        self.__filepath = filepath
        self.__default_resolution = resolutions[0]
        self.__all_resolutions = resolutions
        self.__norm = norm
        self.chrom_dot_sizes = strawC.getChromosomes(filepath)
        self.__footer = {}
        for chromosome in self.chrom_dot_sizes:
            chrom = chromosome.name
            if chrom.lower() == 'all':
                continue
            print('Getting norm vectors for', chrom, flush=True)
            self.__footer[chrom] = {}
            for res in resolutions:
                self.__footer[chrom][res] = strawC.getNormExpVectors(filepath, chrom, chrom, "observed",
                                                                     norm, "BP", res)

    def grab_intra_records(self, chrom1: str, cx1: int, cx2: int, cy1: int, cy2: int, resolution: int):
        footer = self.__footer[chrom1][resolution]
        return strawC.getRecords(self.__filepath, cx1, cx2, cy1, cy2,
                                 resolution, footer.foundFooter, footer.version,
                                 footer.c1, footer.c2, footer.numBins1, footer.numBins2,
                                 footer.myFilePos, footer.unit, footer.norm,
                                 footer.matrixType, footer.c1Norm, footer.c2Norm,
                                 footer.expectedValues)

    def grab_intra_region(self, chrom1: str, cx1: int, cx2: int, cy1: int, cy2: int, resolution: int, width: int):
        result = self.grab_intra_records(chrom1, cx1, cx2, cy1, cy2, resolution)
        row_indices, col_indices, data = list(), list(), list()
        for record in result:
            if cx1 <= record.binX <= cx2 and cy1 <= record.binY <= cy2:
                row_indices.append(record.binX)
                col_indices.append(record.binY)
                data.append(record.counts)
            if record.binX != record.binY and cx1 <= record.binY <= cx2 and cy1 <= record.binX <= cy2:
                row_indices.append(record.binY)
                col_indices.append(record.binX)
                data.append(record.counts)

        row_indices = (np.asarray(row_indices) - cx1) / resolution
        col_indices = (np.asarray(col_indices) - cy1) / resolution
        matrix = sparse.coo_matrix((data, (row_indices.astype(int), col_indices.astype(int))),
                                   shape=(width, width)).toarray()
        matrix[np.isnan(matrix)] = 0
        matrix[np.isinf(matrix)] = 0
        return matrix

    def get_data_from_coordinates(self, coordinates, width: int, resolution: int):
        (chrom1, x1, y1) = coordinates
        cx1 = x1 * resolution
        cx2 = (x1 + width - 1) * resolution
        cy1 = y1 * resolution
        cy2 = (y1 + width - 1) * resolution
        return self.grab_intra_region(chrom1, cx1, cx2, cy1, cy2, resolution, width)
