/// These are the selection schemes we are using for this project
/// These selection functions are selecting on the assumption that the problem is a maximization problem

#ifndef SEL_H
#define SEL_H

///< standard headers
#include <algorithm>
#include <map>
#include <utility>
#include <cmath>
#include <numeric>
#include <set>

///< empirical headers
#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

class Selection
{
  // object types we are using in this class
  public:
    // vector of any ids
    using ids_t = emp::vector<size_t>;
    // vector type of org score
    using phenotype_t = emp::vector<double>;
    // matrix type of org with multiple scores
    using fmatrix_t = emp::vector<phenotype_t>;
    // vector holding population genomes
    using gmatrix_t = emp::vector<emp::vector<double>>;
    // map holding population id groupings by fitness (keys in decending order)
    using fitgp_t = std::map<double, ids_t, std::greater<double>>;
    // sorted score vector w/ position id and score
    using sorted_t = emp::vector<std::pair<size_t,double>>;
    // vector of double vectors for K neighborhoods
    using neigh_t = emp::vector<phenotype_t>;
    // vector of vector size_t for Pareto grouping
    using pareto_g_t = emp::vector<emp::vector<size_t>>;
    // pairing of position_id and fitness associated with id
    using id_f_t = std::pair<size_t,double>;


  public:

    Selection(emp::Ptr<emp::Random> rng = nullptr) : random(rng) {emp_assert(rng);}

    ///< sanity check functions

    // make sure all ids are unique
    template <class T>
    bool ParetoUnique(const emp::vector<emp::vector<T>> & v, const size_t & size)
    {
      std::set<T> check;
      for(const emp::vector<T> & u : v)
      {
        for(T x : u)
        {
          check.insert(x);
        }
      }

      return size == check.size();
    }

    bool ParetoNonZero(const phenotype_t & score)
    {
      for(const auto & s : score)
      {
        if(s==0.0)
        {
          return false;
        }
      }
      return true;
    }

    ///< helper functions

    // distance function between two values
    double Distance(double a, double b) {return std::abs(a - b);}

    // print vectors
    template <class T>
    void PrintVec(const emp::vector<T> &v, const std::string s)
    {
      std::cout << s << ": ";
      for(auto x : v) {std::cout << (double) x << ",";}
      std::cout << std::endl;
    }


    ///< population grouping

    /**
     * Fitness Group Structure:
     *
     * This function will return a map of double,vector pairing where
     * each vector corresponds to solution ids share the same fitness
     *
     * @param score Vector containing all solution scores.
     *
     * @return Map grouping population ids by fitness (keys (fitness) in decending order)
     */
    fitgp_t FitnessGroup(const phenotype_t & score);


    ///< selector functions

    /**
     * Truncation Selector:
     *
     * This function finds the top 't' performing solutions.
     * Once the top 't' solutions are found, each will be placed pop size / t times in the parent list.
     *
     * @param top Number of top performing solutions to pick.
     * @param N Population size.
     * @param group Map with fitness,vector id pairs (descending order expected). Doesn't have to be fitness but some double value in best to worst order.
     *
     * @return Vector with parent id's that are selected.
     */
    ids_t Truncation(const size_t top, const size_t N, const fitgp_t & group);

    /**
     * Tournament Selector:
     *
     * This function holds tournaments for 't' solutions randomly picked from the population.
     * In the event of ties, a solution will be selected randomly from all solutions that tie.
     *
     * @param t Tournament size.
     * @param score Vector holding the population score.
     *
     * @return Parent id of solution that won the tournament.
     */
    size_t Tournament(const size_t t, const phenotype_t & score);

    /**
     * Epsilon Lexicase Selector:
     *
     * This function will iterate through individual testcases.
     * While filtering down the population on each one, only taking the top performers.
     * The top performers must be some within some distance of the top performance to pass.
     *
     * In the event of ties on the last testcase being used, a solution will be selected at random.
     *
     * @param mscore Matrix of solution fitnesses (must be the same amount of fitnesses per solution)(mscore.size => # of orgs).
     * @param epsi Epsilon threshold value.
     * @param M Number of traits we are expecting.
     *
     * @return A single winning solution id.
     */
    size_t EpsiLexicase(const fmatrix_t & mscore, const double epsi, const size_t M);

  private:

    // random pointer from world.h
    emp::Ptr<emp::Random> random;
};

///< population groupings

Selection::fitgp_t Selection::FitnessGroup(const phenotype_t & score)
{
  // quick checks
  emp_assert(0 < score.size());

  // place all solutions in map based on score
  fitgp_t group;
  for(size_t i = 0; i < score.size(); ++i)
  {
    // didn't find in group
    if(group.find(score[i]) == group.end())
    {
      ids_t p{i};
      group[score[i]] = p;
    }
    else{group[score[i]].push_back(i);}
  }

  return group;
}



///< selector functions

Selection::ids_t Selection::Truncation(const size_t top, const size_t N, const fitgp_t & group)
{
  // quick checks
  emp_assert(0 < top); emp_assert(0 < N);
  emp_assert(top <= N); emp_assert(0 < group.size());

  // in the event that we are asking for the whole population
  // just return the vector of ids from 0 to pop size (N)
  if(top == N)
  {
    ids_t pop;
    for(const auto & g : group)
    {
      for(const auto & id : g.second)
      {
        pop.push_back(id);
      }
    }

    emp_assert(pop.size() == N);
    return pop;
  }

  // go through the ordered scores and get our top solutions
  ids_t topmu;
  for(auto & g : group)
  {
    auto gt = g.second;
    emp::Shuffle(*random, gt);
    for(auto id : gt)
    {
      topmu.push_back(id);
      if(topmu.size() == top) {break;}
    }

    if(topmu.size() == top) {break;}
  }

  // insert the correct amount of ids
  ids_t parent;
  size_t lm = N / top;
  for(auto id : topmu)
  {
    for(size_t i = 0; i < lm; ++i){parent.push_back(id);}
  }

  emp_assert(parent.size() == N);

  return parent;
}

size_t Selection::Tournament(const size_t t, const phenotype_t & score)
{
  // quick checks
  emp_assert(0 < t); emp_assert(0 < score.size());
  emp_assert(t <= score.size());

  // get tournament ids
  emp::vector<size_t> tour = emp::Choose(*random, score.size(), t);

  // store all scores for the tournament
  emp::vector<double> subscore(t);
  for(size_t i = 0; i < tour.size(); ++i) {subscore[i] = score[tour[i]];}

  // group scores by fitness and position;
  auto group = FitnessGroup(subscore);
  emp::vector<size_t> opt = group.begin()->second;

  //shuffle the vector with best fitness ids
  emp::Shuffle(*random, opt);
  emp_assert(0 < opt.size());

  return tour[opt[0]];
}

size_t Selection::EpsiLexicase(const fmatrix_t & mscore, const double epsi, const size_t M)
{
  // quick checks
  emp_assert(0 < mscore.size()); emp_assert(0 <= epsi); emp_assert(0 < M);

  // create vector of shuffled testcase ids
  ids_t test_id(M);
  std::iota(test_id.begin(), test_id.end(), 0);
  emp::Shuffle(*random, test_id);

  // vector to hold filterd elite solutions
  ids_t filter(mscore.size());
  std::iota(filter.begin(), filter.end(), 0);

  // iterate through testcases until we run out or have a single winner
  size_t tcnt = 0;
  while(tcnt < M && filter.size() != 1)
  {
    // testcase we are randomly evaluating
    size_t testcase = test_id[tcnt];

    // create vector of current filter solutions
    phenotype_t scores(filter.size());
    for(size_t i = 0; i < filter.size(); ++i)
    {
      // make sure each solutions vector is the right size
      emp_assert(mscore[filter[i]].size() == M);
      scores[i] = mscore[filter[i]][testcase];
    }

    // group org ids by performance in descending order
    fitgp_t group = FitnessGroup(scores);

    // update the filter vector with pop ids that are worthy
    ids_t temp_filter = filter;
    filter.clear();
    for(const auto & p : group)
    {
      if(Distance(group.begin()->first, p.first) <= epsi)
      {
        for(auto id : p.second)
        {
          filter.push_back(temp_filter[id]);
        }
      }
      else{break;}
    }

    ++tcnt;
  }

  // Get a random position from the remaining filtered solutions (may be one left too)
  emp_assert(0 < filter.size());
  size_t wid = emp::Choose(*random, filter.size(), 1)[0];

  return filter[wid];
}

#endif