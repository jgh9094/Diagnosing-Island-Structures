/// Set of diagnostic problems being used in research project
#ifndef STRUCTURE_H
#define STRUCTURE_H

///< standard headers
#include <algorithm>
#include <utility>

///< empirical headers
#include "emp/base/vector.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

class Structure
{
	public:

		// typename for structure structure[island#][world_position]
		using structure_t = emp::vector<emp::vector<size_t>>;
		// typename for island, each value is a position in the world
		using island_t = emp::vector<size_t>;

	public:

		Structure(const size_t cnt_, const size_t sze_, emp::Ptr<emp::Random> rng = nullptr) : isl_count(cnt_), isl_size(sze_), random(rng)
		{
			// generate number of islands & populate each island with dummy numbers
			structure.resize(isl_count);

			// initialize with position id that is outside of possible range
			for(island_t & isl : structure) {isl.resize(isl_size, cnt_ * sze_);}

			// initialize island position counters
			island_cntrs.resize(cnt_, 0);
		}

		Structure() = delete;
		Structure(const Structure &) = delete;
		Structure(Structure &&) = delete;
		~Structure() { ; }
		Structure &operator=(const Structure &) = default;
		Structure &operator=(Structure &&) = default;

		///< Getters
		const structure_t & GetStructure() const {return structure;}
		const island_t & GetIslandCounters() const {return island_cntrs;}
		const island_t & GetIsland(const size_t & isl) const {emp_assert(0 <= isl); emp_assert(isl < isl_count); return structure[isl];}
		const size_t & GetIslandID(const size_t & isl, const size_t & pos) const
		{
			// quick checks
			emp_assert(0 <= isl); emp_assert(isl < isl_count);
			emp_assert(0 <= pos); emp_assert(pos < isl_size);

			return structure[isl][pos];
		}

		const size_t & GetIslandCount() const {return isl_count;}
		const size_t & GetIslandSize() const {return isl_size;}

		///< Helpers
		void PrintIslands();

		/**
		 * @brief Migration will occur and swap solutions between islands
		 *
		 * @param world each the island each solution in the world belongs too
		 */
		void GenerateStructure(const island_t & population);

		/**
		 * @brief Will call migration protocol for specified island count
		 */
		void Migration(const size_t & num);

		/**
		 * @brief Will call migration protocol for four islands
		 */
		void MigrationFourIsls(const size_t & num);

		/**
		 * @brief Guard to check that each island has the expected number of solutions
		 *
		 * @return bool if passes or not
		 */
		bool GRD_IslandSze();

		/**
		 * @brief Guard to check that there are correct number of islands
		 *
		 * @return bool if passes or not
		 */
		bool GRD_IslandCnt();

		/**
		 * @brief Guard to check that each solution's associated island is within range
		 *
		 * @return bool if passes or not
		 */
		bool GRD_IslandRange();

		/**
		 * @brief Guard to check that non solution is overlapping between islands
		 *
		 * @return bool if passes or not
		 */
		bool GRD_IslandUnique();

		/**
		 * @brief Guard to check that each solution in an island is accounted for
		 *
		 * @return bool if passes or not
		 */
		bool GRD_EveryoneAccounted();

		/**
		 * @brief Guard to check lower bound of of population island ids
		 *
		 * @return bool if passes or not
		 */
		bool GRD_PopulationIdsLower(const island_t & population) {return *std::min_element(population.begin(), population.end()) == 0;}

		/**
		 * @brief Guard to check upper bound of of population island ids
		 *
		 * @return bool if passes or not
		 */
		bool GRD_PopulationIdsUpper(const island_t & population) {return *std::max_element(population.begin(), population.end())  == isl_count - 1;}

	private:

		// number of islands
		size_t isl_count = 0;
		// number of organisms on each island
		size_t isl_size = 0;
		// island strucuture => structure [ISLANDS] [SOLUTION_POPULATION_POS]
		structure_t structure;
        // island insertion counters
        island_t island_cntrs;
		// enums for island counts
		enum IslandCount {TWO=2,FOUR=4,EIGHT=8,SIXTEEN=16,THIRTYTWO=32,SIXTYFOUR=64};
		// random pointer from world.h
		emp::Ptr<emp::Random> random;
};

void Structure::GenerateStructure(const island_t & population)
{
    // checks
    emp_assert(population.size() == isl_count * isl_size);
    // make sure that the last island is in the population
    emp_assert(GRD_PopulationIdsUpper(population));
    // make sure that the first island is in the population
    emp_assert(GRD_PopulationIdsLower(population));

	// reset the island counters
	std::fill(island_cntrs.begin(), island_cntrs.end(), 0);

    // populate the current structure
    for(size_t i = 0; i < population.size(); ++i)
    {
		// island we are adding to
		const size_t isl = population[i];

        // add solution isl_id to appropiate island
        structure[isl][island_cntrs[isl]] = i;

        // increment island counter
        island_cntrs[isl]++;
    }

    // make sure all guards hold;
	emp_assert(GRD_IslandSze()); emp_assert(GRD_IslandCnt()); emp_assert(GRD_IslandRange());
	emp_assert(GRD_IslandUnique()); emp_assert(GRD_EveryoneAccounted());
}

void Structure::Migration(const size_t & num)
{
	// make sure all guards hold;
	emp_assert(GRD_IslandSze()); emp_assert(GRD_IslandCnt()); emp_assert(GRD_IslandRange());
	emp_assert(GRD_IslandUnique()); emp_assert(GRD_EveryoneAccounted());

	switch (isl_count)
	{
		case static_cast<size_t>(IslandCount::FOUR):
			MigrationFourIsls(num);
			break;

		default:
			std::cout << "ERROR UNKNOWN ISLAND COUNT" << std::endl;
			emp_assert(false);
			break;
	}
}

void Structure::MigrationFourIsls(const size_t & num)
{
	// make sure all guards hold
	emp_assert(isl_count == static_cast<size_t>(IslandCount::FOUR)); emp_assert(GRD_IslandSze());
	emp_assert(GRD_IslandCnt()); emp_assert(GRD_IslandRange()); emp_assert(GRD_IslandUnique());
	emp_assert(GRD_EveryoneAccounted());

	// island swaps: 0 <=> 3 & 1 <=> 2 & 0 <=> 1 & 2 <=> 3
	// this ordering of swaps guarentees that no solution can return to the same island
	const emp::vector<std::pair<int, int>> swaps = {{0,3},{1,2},{0,1},{2,3}};

	// perform swaps
	for(const auto & p : swaps)
	{
		// get ids we are swapping between islands
		const emp::vector<size_t> grants = emp::Choose(*random, isl_size, num);

		// swap solutions between islands
		for(const size_t & g : grants) {std::swap(structure[p.first][g], structure[p.second][g]);}
	}
}


bool Structure::GRD_IslandSze()
{
	for(const island_t & isl : structure)
	{
		if (isl_size != isl.size()) {return false;}
	}
	return true;
}

bool Structure::GRD_IslandCnt()
{
	return structure.size() == isl_count;
}

bool Structure::GRD_IslandUnique()
{
	emp::vector<bool> uni(isl_count * isl_size, false);
	for(const island_t & isl : structure)
	{
		for(const size_t & id : isl)
		{
			if(uni[id]) {return false;}
			uni[id] = true;
		}
	}
	return true;
}

bool Structure::GRD_IslandRange()
{
	for(const island_t & isl : structure)
	{
		for(const size_t & id : isl)
		{
			if(id < 0 or isl_count * isl_size <= id) {return false;}
		}
	}
	return true;
}

bool Structure::GRD_EveryoneAccounted()
{
	emp::vector<bool> acc(isl_count * isl_size, false);
	for(const island_t & isl : structure)
	{
		for(const size_t & id : isl) {acc[id] = true;}
	}
	return std::accumulate(acc.begin(), acc.end(), 0) == isl_count * isl_size;
}


void Structure::PrintIslands()
{
	for(size_t isl = 0; isl < isl_count; ++isl)
	{
		std::cout << "Island " << isl << ": ";
		for(const size_t & id : structure[isl]) {std::cout << id << ",";}
		std::cout << std::endl;
	}
}

#endif