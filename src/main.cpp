// test code
#include "csv.cpp"
#include "matrix.cpp"
#include <fstream>
#include <iostream>
#include <stdlib.h> // atof
#include <vector>

// reduce base vectors (subtract sub-genotypes from covering super-genotypes)
template<typename T>
void reduce_row_vectors(matrix<T>&);
// read CSV input into matrix in transposed form
template<typename T>
void read_csv_transposed(std::istream&, matrix<T>&);
// Compile a set of base vectors for a finite vector space,
// given a RREF-transformed version of the matrix.
template<typename T>
void extract_base_vectors(const matrix<T>&, const matrix<T>&, std::vector<int>&, matrix<T>&);
template<typename T, typename S>
void print_matrix(const matrix<T>&, const std::vector<S>);

void print_usage(char *prog_name)
{
  std::cerr << "usage: " << prog_name << " haplotypes.csv" << std::endl;
}

int main(int argc, char *argv[])
{
  matrix<double> mat_2d = matrix<double>();

  // check command line params
  if (argc < 2)
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  // read CSV file (transposed)
  std::ifstream file(argv[1]);
  read_csv_transposed(file, mat_2d);

  // make a copy of the input haplotypes, then transform matrix to RREF
  matrix<double> mat_haplotypes = mat_2d;
  to_reduced_row_echelon_form(mat_2d);

  // extract base vectors
  std::vector<int> idx_base;
  matrix<double> mat_base = matrix<double>();
  extract_base_vectors(mat_2d, mat_haplotypes, idx_base, mat_base);
  std::cout << std::endl << "Base contains " << idx_base.size() << " vectors:" << std::endl;
  print_matrix(mat_base, idx_base);

  // reduce vector base to avoid subtractions
  reduce_row_vectors(mat_base);
  std::cout << std::endl << "Reduced base:" << std::endl;
  print_matrix(mat_base, idx_base);

  return EXIT_SUCCESS;
}


/** Reads data from a CSV-formatted file.
 *  Expects first row to contain the header, first column to contain IDs. */
template<typename T>
void read_csv_transposed(std::istream& input, matrix<T>& M)
{
  CSVRow row;
  input >> row; // read header row
  M.rows = row.size()-1; // first column expected to contain IDs
  M.data = std::vector<std::vector<double> >(M.rows);
  while (input >> row)
  {
    for (int i=0; i<M.rows; ++i)
    {
      double cell = std::atof(row[i+1].c_str());
      M.data[i].push_back(cell);
    }
    M.columns++;
  }
}

/** Extract vector base that generates a given vector space. */
template<typename T>
void extract_base_vectors(const matrix<T> &M, const matrix<T> &D, std::vector<int> &idx, matrix<T> &V)
{
  // find pivot columns for each row (-> indices of base vectors)
  for (int i=0; i<M.rows; ++i)
  {
    int j=0;
    while (j<M.columns && M.data[i][j]==0)
      j++;
    bool is_last_one = (j==M.columns-1 && M.data[i][j]==1);
    if (j<M.columns-1 || is_last_one)
      idx.push_back(j);
  }

  // compile vector base (transpose to facilitate row subtractions)
  for (unsigned i=0; i<idx.size(); ++i)
  {
    V.data.push_back(std::vector<double>(M.rows));
    for (int j=0; j<D.rows; ++j)
      V.data[i][j] = D.data[j][idx[i]];
  }
  assert(V.data.size() == idx.size());
  assert((int)(V.data[0].size()) == D.rows);
  V.rows = (int)(idx.size());
  V.columns = D.rows;
}

template<typename T>
void reduce_row_vectors(matrix<T>& M)
{
  bool reiterate = true;
  while (reiterate)
  {
    reiterate = false;
    for (int v=0; v<M.rows; ++v)
      for (int w=0; w<M.rows; ++w)
        if (v!=w)
        {
          bool covered = true;
          for (int i=0; i<M.columns; ++i)
            if (M.data[v][i] > M.data[w][i])
            {
              covered = false;
              break;
            }
          if (covered)
          {
            add_multiple_row(M, w, v, -1);
            reiterate = true;
          }
        }
  }
}

template<typename T, typename S>
void print_matrix(const matrix<T>& M, const std::vector<S> id)
{
  assert((int)(id.size()) == M.rows);
  for (unsigned i=0; i<id.size(); ++i)
  {
    std::cout << id[i] << ":\t";
    for (int j=0; j<M.columns-1; ++j)
      std::cout << M.data[i][j] << ",";
    std::cout << M.data[i][M.columns-1] << std::endl;
  }
  std::cout << std::endl;
}
