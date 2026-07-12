#include "spglib.h"

#include <math.h>

enum mock_case {
  MOCK_UNKNOWN,
  MOCK_IDENTITY,
  MOCK_DIAMOND,
  MOCK_TRANSLATED,
  MOCK_ZERO_COUNT,
  MOCK_API_FAILURE,
  MOCK_INVALID_GROUP,
  MOCK_EXCESSIVE_COUNT,
  MOCK_CAPACITY
};

static int near(double left, double right) {
  return fabs(left - right) < 1.0e-12;
}

static enum mock_case classify(const double lattice[3][3],
                               const double position[][3], const int types[],
                               int num_atom) {
  if (num_atom == 2 && types[0] == 6 && types[1] == 6 &&
      near(position[1][0], 0.25) && near(position[1][1], 0.25) &&
      near(position[1][2], 0.25) && near(lattice[0][0], 0.0) &&
      near(lattice[0][1], 0.5) && near(lattice[0][2], 0.5) &&
      near(lattice[1][0], 0.5)) {
    return MOCK_DIAMOND;
  }
  if (num_atom == 2 && types[0] == 14 && types[1] == 14 &&
      near(position[0][0], 0.5) && near(position[1][0], 0.0) &&
      near(lattice[1][0], 0.4) && near(lattice[1][1], 1.5)) {
    return MOCK_TRANSLATED;
  }
  if (num_atom == 3 && types[0] == 14 && types[1] == 8 && types[2] == 1 &&
      near(position[0][0], 0.11) && near(position[1][1], 0.07) &&
      near(position[2][2], 0.19) && near(lattice[1][0], 0.2) &&
      near(lattice[2][0], 0.1) && near(lattice[2][1], 0.3)) {
    return MOCK_IDENTITY;
  }
  if (num_atom == 1 && types[0] == 4 && near(position[0][0], 0.04)) {
    return MOCK_ZERO_COUNT;
  }
  if (num_atom == 1 && types[0] == 5 && near(position[0][0], 0.05)) {
    return MOCK_API_FAILURE;
  }
  if (num_atom == 1 && types[0] == 6 && near(position[0][0], 0.0)) {
    return MOCK_INVALID_GROUP;
  }
  if (num_atom == 1 && types[0] == 7 && near(position[0][0], 0.07)) {
    return MOCK_EXCESSIVE_COUNT;
  }
  if (num_atom == 1 && types[0] == 8) return MOCK_CAPACITY;
  return MOCK_UNKNOWN;
}

static void set_identity(int rotation[3][3], double translation[3]) {
  int row;
  int column;

  for (row = 0; row < 3; ++row) {
    translation[row] = 0.0;
    for (column = 0; column < 3; ++column) {
      rotation[row][column] = row == column ? 1 : 0;
    }
  }
}

int spg_get_multiplicity(const double lattice[3][3],
                         const double position[][3], const int types[],
                         int num_atom, double symprec) {
  (void)symprec;
  switch (classify(lattice, position, types, num_atom)) {
    case MOCK_UNKNOWN:
      return 0;
    case MOCK_IDENTITY:
      return 1;
    case MOCK_ZERO_COUNT:
      return 0;
    case MOCK_EXCESSIVE_COUNT:
      return 4097;
    default:
      return 2;
  }
}

int spg_get_symmetry(int rotation[][3][3], double translation[][3],
                     int max_size, const double lattice[3][3],
                     const double position[][3], const int types[],
                     int num_atom, double symprec) {
  enum mock_case test_case = classify(lattice, position, types, num_atom);
  int count = spg_get_multiplicity(lattice, position, types, num_atom, symprec);

  if (test_case == MOCK_API_FAILURE || max_size < count) return 0;
  set_identity(rotation[0], translation[0]);
  if (count == 1) return 1;

  set_identity(rotation[1], translation[1]);
  if (test_case == MOCK_DIAMOND) {
    rotation[1][0][0] = -1;
    rotation[1][1][1] = -1;
    rotation[1][2][2] = -1;
    translation[1][0] = 0.25;
    translation[1][1] = 0.25;
    translation[1][2] = 0.25;
  } else if (test_case == MOCK_TRANSLATED) {
    translation[1][0] = 0.5;
  } else if (test_case == MOCK_INVALID_GROUP || test_case == MOCK_CAPACITY) {
    rotation[1][0][0] = 0;
    rotation[1][0][1] = -1;
    rotation[1][1][0] = 1;
    rotation[1][1][1] = 0;
  }
  return count;
}
