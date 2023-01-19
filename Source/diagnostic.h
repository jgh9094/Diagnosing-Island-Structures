/// Set of diagnostic problems being used in research project
#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

///< standard headers
#include <algorithm>

///< empirical headers
#include "emp/base/vector.hpp"

class Diagnostic
{
    public:
        //typename for target vector
        using phenotype_t = emp::vector<double>;
        using genome_t = emp::vector<double>;
        using optimal_t = emp::vector<bool>;

  public:

        Diagnostic(double t_, double b_) : target(t_), lower_bnd(b_) {lbnd_set = true;}

        Diagnostic() {;}

        // getters
        double GetTarget() const {return target;}
        double GetCredit() const {return lower_bnd;}

        //setters
        void SetTarget(double t_) {target = t_;}
        void SetCredit(double b_) {lbnd_set = true; lower_bnd = b_;}


        ///< Functions that deal with diagnostic problem scoring

        /**
         * Exploitation function:
         *
         * All genes in genome evaluated treated as independent optimization tasks.
         * In other words, there is no interactions between genes.
         * The phenotype vector is just the genome, with the caveat that genes cannot go over associated target value.
         *
         * @param g genome from organism being evaluated.
         *
         * @return phenotype vector that is calculated from 'g'.
         */
        phenotype_t ExploitationRate(const genome_t & g);

        /**
         * Ordered Exploitation function:
         *
         * All genes in genome expected to follow descending order for evaluation.
         * While in decending order, phenotype[i] = genome[i].
         * If descending order is broken, max_credit is assigned to each following trait .
         * Note that this problem explicitly forces optimization from start to end.
         *
         * @param g genome from organism being evaluated.
         *
         * @return phenotype vector that is calculated from 'g'.
         */
        phenotype_t OrderedExploitation(const genome_t & g);

        /**
         * Contradictory Objectives function:
         *
         * Solutions are pressured to optimize one gene, while all other genes are ignored and set to max_credit.
         * To calculate phenotype vector, first we find the maximum gene value at position i.
         * After, phenotype[i] = genome[i], and every every other trait is set to max_credit.
         *
         * @param g genome from organism being evaluated.
         *
         * @return phenotype vector that is calculated from 'g'.
         */
        phenotype_t ContradictoryObjectives(const genome_t & g);

        /**
         * Multi-path Exploration function:
         *
         * Solutions are must be able to simultaneously explore multiple pathways, with the goal of pursuing the one that leads to the optimum.
         * To calculate phenotype vector, first we find the maximum gene value at position i.
         * From that gene forward, while genes follow decending order, phenotype[i] = genome[i].
         * This process stops until the decending order is broken, or the end of the genome is reached.
         * Note, the optimum is reached when i = 0, the start of the genome.
         *
         * @param g genome from organism being evaluated.
         * @return phenotype vector that is calculated from 'g'.
         */
        phenotype_t MultiPathExploration(const genome_t & g);

        /**
         * Multi-valley Crossing function:
         *
         * Solutions are pressured to cross valleys of different widths at each gene.
         * We use a peaks vector that supplied to calculate the penalty at each valley.
         * We use the difference between the gene value and penalty value to determine
         * the reduction when calculating the trait value.
         *
         * phenotype[i] = peaks[ floor(genome[i]) ] -  (genome[i] - peaks[ floor(genome[i]) ])
         *
         * @param g genome from organism being evaluated.
         * @param peaks vector of peaks for each floored gene value
         * @param dips_start gene value where peaks begin
         * @param dips_end gene value where dips end
         *
         * @return phenotype vector that is calculated from 'g'.
         */
        phenotype_t MultiValleyCrossing(const genome_t & g, const phenotype_t & peaks, const double & dips_start, const double & dips_end);

        ///< Functions that deal with interpretation of phenotype vectors

        /**
         * Optimized Vector:
         *
         * Will return a boolean vector that lists what genes are optimized relative to the target vector.
         * A gene is optimized if it meets the target value accuracy % requirement.
         *
         * @param g genome from organism being evaluated.
         * @param acc This value is the accuracy % needed to be considered optimized
         *
         * @return boolean vector that displays optimized traits.
        */
        optimal_t SatisfactoryVector(const genome_t & g, const double & acc);

    private:
    // holds vector of target objective values
    double target = 0.0;
    // holds maximum credit allowed for error
    double lower_bnd = 0.0;
    // lower bound set?
    bool lbnd_set = false;
};

///< diagnostic problem implementations

Diagnostic::phenotype_t Diagnostic::ExploitationRate(const genome_t & g)
{
  // quick checks
  emp_assert(0 < g.size());

  // intialize vector with size g
  phenotype_t phenotype = g;

  return phenotype;
}

Diagnostic::phenotype_t Diagnostic::OrderedExploitation(const genome_t & g)
{
  // quick checks
  emp_assert(g.size() > 0); emp_assert(lbnd_set);

  // find where descending order breaks
  const auto it = std::is_sorted_until(g.begin(), g.end(), std::greater<>());

  // if sorted, return same vector
  if(it == g.end())
  {
    phenotype_t phenotype = g;
    return phenotype;
  }
  // else fill in appropiately
  else
  {
    // initialize phenotype to the same size
    phenotype_t phenotype(g.size());

    // calculate cutoff point where descending order is broken
    const size_t cutoff = std::distance(g.begin(), it);

    // everything up to unsorted
    for(size_t i = 0; i < cutoff; ++i) {phenotype[i] = g[i];}

    // everything after unsorted
    for(size_t i = cutoff; i < phenotype.size(); ++i) {phenotype[i] = lower_bnd;}

    return phenotype;
  }
}

Diagnostic::phenotype_t Diagnostic::ContradictoryObjectives(const genome_t & g)
{
  // quick checks
  emp_assert(g.size() > 0);

  // intialize phenotype vector
  phenotype_t phenotype(g.size());

  // find max value position
  const size_t max_gene = std::distance(g.begin(), std::max_element(g.begin(), g.end()));

  // set all phenotype vector values
  for(size_t i = 0; i < phenotype.size(); ++i)
  {
    if(i == max_gene) {phenotype[i] = g[max_gene];}
    else {phenotype[i] = lower_bnd;}
  }

  return phenotype;
}

Diagnostic::phenotype_t Diagnostic::MultiPathExploration(const genome_t & g)
{
  // quick checks
  emp_assert(0 < g.size()); emp_assert(lbnd_set);

  // intialize vector with size g
  phenotype_t phenotype(g.size());

  // find max value position
  const auto opti_it = std::max_element(g.begin(), g.end());
  const size_t opti = std::distance(g.begin(), opti_it);

  // find where order breaks
  const size_t sort = std::distance(g.begin(), std::is_sorted_until(opti_it, g.end(), std::greater<>()));

  // left of optimal value found
  for(size_t i = 0; i < opti; ++i) {phenotype[i] = lower_bnd;}
  // middle of optimal value till order broken
  for(size_t i = opti; i < sort; ++i) {phenotype[i] = g[i];}
  // right of order broken
  for(size_t i = sort; i < phenotype.size(); ++i) {phenotype[i] = lower_bnd;}

  return phenotype;
}

Diagnostic::phenotype_t Diagnostic::MultiValleyCrossing(const phenotype_t & p, const phenotype_t & peaks, const double & dips_start, const double & dips_end)
{
  // quick checks
  emp_assert(p.size() > 0); emp_assert(lbnd_set);
  emp_assert(peaks.size() > 0);

  // intialize vector with size p
  phenotype_t phenotype(p.size());

  for(size_t i = 0; i < p.size(); ++i)
  {
    if (p[i] <= dips_start || p[i] >= dips_end)
    {
      phenotype[i] = p[i];
    }
    else
    {
      phenotype[i] = 2.0 * peaks[static_cast<size_t>(p[i])] - p[i];
    }
  }

  return phenotype;
}

///< phenotype vector interpretation implementations

Diagnostic::optimal_t Diagnostic::SatisfactoryVector(const genome_t & g, const double & acc)
{
  // quick checks
  emp_assert(g.size() > 0);
  emp_assert(0.0 < target);
  emp_assert(0.0 < acc);
  emp_assert(acc <= 1.0);

  // initialize optimized vector with all false
  optimal_t optimize(g.size(), false);
  const double threshold = acc * target;

  // iterate through genome and check optimality
  for(size_t i = 0; i < g.size(); ++i)
  {
    if(threshold <= g[i]) {optimize[i] = true;}
  }

  return optimize;
}

#endif