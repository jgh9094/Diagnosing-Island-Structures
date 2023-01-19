/// organims class that holds its genome and phenotype, both numerical vectors of a predetermined dimensionality

#ifndef ORG_H
#define ORG_H

///< standard headers
#include <algorithm>

///< empirical headers
#include "emp/base/vector.hpp"

///< coordiante we start from
constexpr double START_LB = 0.0;

class Org
{
	public:
		// genome vector type
		using genome_t = emp::vector<double>;
		// phenotype vector type
		using phenotype_t = emp::vector<double>;
		// optimal gene vector type
		using optimal_t = emp::vector<bool>;

	public:

		// every org after starting generation
		Org(genome_t _g)
		{
			// quick checks
			emp_assert(genome.size() == 0); emp_assert(M == 0); emp_assert(!evaluated);

			// set all requried variables
			M = _g.size();
			start_pos = _g.size();
			streak = _g.size();
			genome = _g;
		}

		// ask Charles
		Org(const Org &) = delete;
		Org(Org &&) = delete;
		~Org() { ; }
		Org &operator=(const Org &) = default;
		Org &operator=(Org &&) = default;

		///< getters

		const genome_t & GetGenome() const {emp_assert(0 < genome.size()); return genome;}
		genome_t & GetGenome() {emp_assert(0 < genome.size()); return genome;}

		const phenotype_t & GetPhenotype() const {emp_assert(evaluated); return phenotype;}
		phenotype_t & GetPhenotype() {emp_assert(evaluated); return phenotype;}

		const optimal_t & GetSatisfactoryVec() const {emp_assert(opti); return optimal;}
		optimal_t & GetSatisfactoryVec() {emp_assert(opti); return optimal;}

		const double & GetAggregate() const {emp_assert(aggregated); return aggregate;}
		double & GetAggregate() {emp_assert(aggregated); return aggregate;}

		// get are we a clone bool
		const bool & GetClone() const {emp_assert(0 < genome.size()); return clone;}
		// get satisfactory trait count
		const size_t & GetCount() const {emp_assert(counted); return count;}
		// get gene count
		const size_t GetM() {emp_assert(0 < M); return M;}
		// get start position
		const size_t & GetStart() const {emp_assert(start_pos != M); return start_pos;}
		// get streak count
		const size_t & GetStreak() const {emp_assert(streak != M); return streak;}
		// Are we optimized at this trait?
		const bool OptimizedAt(const size_t obj) const;
		// get island id
		const size_t GetIslandID() const {emp_assert(islanded); return island;}

		// get evaluated bool
		const bool & GetEvaluated() const {return evaluated;}
		// get optimal bool
		const bool & GetOpti() const {return opti;}
		// get aggregated bool
		const bool & GetAggregated() const {return aggregated;}
		// get counted bool
		const bool & GetCounted() const {return counted;}
		// get streak bool
		const bool & GetStreaked() const {return streaked;}
		// get max trait
		const double & GetMaxTrait() const
		{
			// quick checks
			emp_assert(start); emp_assert(0 < M);
			emp_assert(phenotype.size() == M);

			return phenotype[start_pos];
		}

		///< setters

		// set phenotype vector (recieved from problem.h in world.h or inherited from parent)
		void SetPhenotype(const phenotype_t & p_)
		{
			// make sure that phenotype vector hasn't been set before.
			emp_assert(!evaluated); emp_assert(p_.size() == M); emp_assert(phenotype.size() == 0); emp_assert(0 < M);
			evaluated = true;
			phenotype = p_;
		}

		// set the optimal gene vector (recieved from problem.h in world.h or inherited from parent)
		void SetOptimal(const optimal_t & o_)
		{
			// make sure that optimal gene vector hasn't been set before.
			emp_assert(!opti); emp_assert(o_.size() == M); emp_assert(optimal.size() == 0); emp_assert(0 < M);
			opti = true;
			optimal = o_;
		}

		// set the optimal gene count (called from world.h or inherited from parent)
		void SetCount(size_t c_)
		{
			emp_assert(!counted); emp_assert(0 < M);
			counted = true;
			count = c_;
		}

		// set the aggregated phenotype (called from world.h or inherited from parent)
		void SetAggregate(double a_)
		{
			emp_assert(!aggregated); emp_assert(0 < M);
			aggregated = true;
			aggregate = a_;
		}

		// set the starting position
		void SetStart(size_t s_)
		{
			emp_assert(!start); emp_assert(0 < M);
			start = true;
			start_pos = s_;
		}

		// set the starting position
		void SetStreak(size_t s_)
		{
			emp_assert(!streaked); emp_assert(0 < M);
			streaked = true;
			streak = s_;
		}

		// set the island id
		void SetIslandID(size_t id_)
		{
			emp_assert(!islanded or island != 512);
			islanded = true;
			island = id_;
		}


		///< functions to calculate scores and related data

		/**
		 * Aggregate Score function:
		 *
		 * First we check to see if the phenotype has been calculated.
		 * If yes, aggregated=True and throw expection.
		 * Else, aggregated=False and we calculate it and return aggregate.
		 *
		 * @return aggregate
		*/
		double AggregateScore();

		/**
		 * Count Optimized function:
		 *
		 * First we check to see if the count has been calculated.
		 * If yes, counted=True and throw expection.
		 * Else, countedd=False and we calculate it and return count.
		 *
		 * @return count
		*/
		size_t CountOptimized();

		/**
		 * Find Starting Position
		 *
		 * Find the staring position in the phenotype vector.
		 */
		size_t StartPosition();


		/**
		 * Find Maximum Streak
		 *
		 * Find the biggest streak.
		 */
		size_t CalcStreak();


		///< functions related to the birth of an organism

		/**
		 * Reset function:
		 *
		 * Will reset all variables in an organims when birth occurs.
		 * Function executes if offspring is not a clone
		*/
		void Reset();

		/**
		 * Inherit function:
		 *
		 * Will pass all info from parent to offspring solution.
		 * Function executes if offpsring is an clone
		 *
		 * @param s phenotype vector recived
		 * @param o optimal gene vector recieved
		 * @param c optimal gene count recieved
		 * @param a aggregate phenotype recieved
		 * @param st start position recieved
		 * @param sr streak count recieved
		 * @param i island id recived
		 *
		*/
		void Inherit(const phenotype_t & s, const optimal_t & o, const size_t c, const double a, const size_t st, const size_t sr, const size_t i);

		/**
		 * Me Clone function:
		 *
		 * Will set the clone variable to true.
		*/
		void MeClone() {emp_assert(0 < M); emp_assert(!clone); clone = true;}

	private:
		// organism genome vector
		genome_t genome;

		// organism phenotype vector
		phenotype_t phenotype;
		// phenotype vector set?
		bool evaluated = false;

		// organims gene optimal vector
		optimal_t optimal;
		// gene optimal vector calculated?
		bool opti = false;

		// optimal gene count
		size_t count = 0;
		// gene optimal vector counted?
		bool counted = false;

		// aggregate phenotype
		double aggregate = 0.0;
		// aggregate calculate?
		bool aggregated = false;

		// streak count
		size_t streak = 0;
		// streak calculated?
		bool streaked = false;

		// Number of genes in genome
		size_t M = 0;

		// starting position
		size_t start_pos;
		// starting position located?
		bool start = false;

		// Are we a clone?
		bool clone = false;

		// island id
		size_t island = 512;
		bool islanded = false;
};

///< getters with extra
const bool Org::OptimizedAt(const size_t obj) const
{
  // quick checks
  emp_assert(0 <= obj); emp_assert(obj < M);
  emp_assert(0 < optimal.size()); emp_assert(M == optimal.size());

  return optimal[obj];
}

///< functions to calculate scores and related data

double Org::AggregateScore()
{
  //quick checks
  emp_assert(!aggregated); emp_assert(0 < M);
  emp_assert(phenotype.size() == M, phenotype.size());

  // calculate the aggregate phenotype and set it
  SetAggregate(std::accumulate(phenotype.begin(), phenotype.end(), 0.0));

  return aggregate;
}

size_t Org::CountOptimized()
{
  //quick checks
  emp_assert(!counted); emp_assert(0 < M); emp_assert(opti);
  emp_assert(optimal.size() == M, optimal.size());

  // calculate total optimal genes and set it
  SetCount(std::accumulate(optimal.begin(), optimal.end(), 0));

  return count;
}

size_t Org::StartPosition()
{
  // quick checks
  emp_assert(!start); emp_assert(0 < M);
  emp_assert(phenotype.size() == M);

  // find max value position
  const auto opti_it = std::max_element(phenotype.begin(), phenotype.end());
  SetStart(std::distance(phenotype.begin(), opti_it));

  return start_pos;
}

size_t Org::CalcStreak()
{
  // quick checks
  emp_assert(!streaked); emp_assert(0 < M);
  emp_assert(phenotype.size() == M);

  // get longest streak
  size_t count = 0;
  size_t max_cnt = 0;
  for(auto & s : phenotype)
  {
    if(s > 0.0)
    {
      count++;
    }
    else
    {
      if(count > max_cnt)
      {
        max_cnt = count;
      }
      count = 0;
    }
  }

  streak = max_cnt;

  return streak;
}

///< functions related to the birth of an organism

void Org::Reset()
{
	// quick checks
	emp_assert(0 < M); emp_assert(0 < genome.size());

	// reset phenotype vector stuff
	phenotype.clear();
	evaluated = false;

	// reset optimal gene vector stuff
	optimal.clear();
	opti = false;

	// reset optimal gene count stuff
	count = 0;
	counted = false;

	// reset aggregate phenotype stuff
	aggregate = 0.0;
	aggregated = false;

	// reset starting position info
	start_pos = genome.size();
	start = false;

	// reset clone var
	clone = false;

	// island id
	island = 512;
	islanded = false;
}

void Org::Inherit(const phenotype_t & s, const optimal_t & o, const size_t c, const double a, const size_t st, const size_t sr, const size_t i)
{
  // quick checks
  emp_assert(0 < M); emp_assert(0 < genome.size()); emp_assert(clone);

  // copy everything into offspring solution
  SetPhenotype(s);
  SetOptimal(o);
  SetCount(c);
  SetAggregate(a);
  SetStart(st);
  SetStreak(sr);
  SetIslandID(i);
}

#endif