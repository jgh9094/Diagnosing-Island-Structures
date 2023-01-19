#ifndef ISL_CONFIG_H
#define ISL_CONFIG_H

// empirical  headers
#include "emp/config/config.hpp"

EMP_BUILD_CONFIG(DiaConfig,
  GROUP(WORLD, "World be setup details."),
  VALUE(POP_SIZE,     size_t,      512,    "Population size."),
  VALUE(MAX_GENS,     size_t,    50000,    "Maximum number of generations."),
  VALUE(SEED,            int,        1,    "Random number seed."),

  GROUP(ISLAND, "Island configurations"),
  VALUE(ISL_CNT,           size_t,       4,    "How many islands are there."),
  VALUE(ISL_SIZE,          size_t,     128,    "Size of each island."),
  VALUE(MIGRATION,           bool,       1,    "Is there migration between islands."),
  VALUE(MIGRATION_INT,     size_t,     500,    "How many generations until a migration occurs."),
  VALUE(MIGRATION_SZE,     size_t,       8,    "How many organims are migrating (swapping between islands)."),

  GROUP(DIAGNOSTICS, "How are the diagnostics setup?"),
  VALUE(LOWER_BND,           double,       0.0,      "Lower bound for random starts."),
  VALUE(UPPER_BND,           double,       1.0,      "Upper bound for random starts."),
  VALUE(TARGET,              double,     100.0,      "Target values that genes must reach."),
  VALUE(ACCURACY,            double,      0.99,      "Accuracy percentage needed to be considered an optimal trait."),
  VALUE(DIMENSIONALITY,      size_t,       100,      "Diagnositc dimensionality."),
  VALUE(SELECTION,           size_t,         0,      "Which selection are we doing? \n0: Truncation\n1: Tournament\n2: Espilon Lexicase"),
  VALUE(DIAGNOSTIC,          size_t,         0,      "Which diagnostic are we doing? \n0: Exploitation Rate\n1: Ordered Exploitation\n"
                                                     "2: Contradictory Objectives \n3: Multi-path Exploration"),
  VALUE(VALLEY_CROSSING,       bool,     false,      "Do we add multi-valley crossing layer to the diagnostics?"),

  GROUP(MUTATIONS, "Mutation details."),
  VALUE(MUTATE_PER,       double,     0.007,        "Probability of genes recieving a mutation."),
  VALUE(MEAN,             double,     0.0,          "Mean of nurmal distribution for point mutations on genes."),
  VALUE(STD,              double,     1.0,          "Standard deviation of normal distribution for point mutations on genes."),

  GROUP(TRUNCATION, "Parameters for truncation selection."),
  VALUE(TRUNC_SIZE,       size_t,     8,            "Truncation size."),

  GROUP(TOURNAMENT, "Parameters for tournament selection."),
  VALUE(TOUR_SIZE,        size_t,     8,            "Tournament size."),

  GROUP(LEXICASE, "Parameters for lexicase selection."),
  VALUE(LEX_EPS,          double,            0.0,       "Epsillon value for differences between best solution when filtering."),

  GROUP(OUTPUT, "Output configurations"),
  VALUE(OUTPUT_DIR,     std::string,              "./",          "Directory where we put all data in")
)

#endif