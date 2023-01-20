/// World that will manage solutions during the evolutionary run

#ifndef DIA_WORLD_H
#define DIA_WORLD_H

///< standard headers
#include <functional>
#include <map>
#include <set>
#include <fstream>
#include <string.h>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>

///< empirical headers
#include "emp/Evolve/World.hpp"
#include "emp/math/random_utils.hpp"

///< experiment headers
#include "config.h"
#include "org.h"
#include "diagnostic.h"
#include "selection.h"
#include "structure.h"

class DiagWorld : public emp::World<Org>
{
  // object types for consistency between working class
  public:
    ///< Org related

    // solution genome + diagnotic problem types
    using genome_t = emp::vector<double>;
    // phenotype vector for a solution
    using phenotype_t = emp::vector<double>;
    // boolean optimal vector per objective
    using optimal_t = emp::vector<bool>;

    ///< selection related types

    // vector of position ids
    using ids_t = emp::vector<size_t>;
    // matrix of population phenotype vectors
    using fmatrix_t = emp::vector<phenotype_t>;
    // matrix of population genomes
    using gmatrix_t = emp::vector<genome_t>;
    // map holding population id groupings by fitness (keys in decending order)
    using fitgp_t = std::map<double, ids_t, std::greater<double>>;

    ///< world related types

    // evaluation function type
    using eval_t = std::function<double(Org &)>;
    // selection function type
    using sele_t = std::function<ids_t()>;

    ///< data tracking stuff (ask about)
    using nodef_t = emp::Ptr<emp::DataMonitor<double>>;
    using nodeo_t = emp::Ptr<emp::DataMonitor<size_t>>;

  public:

    DiagWorld(DiaConfig & _config) : emp::World<Org>("",false), config(_config), data_file(_config.OUTPUT_DIR() + "data.csv")
    {
      emp_assert(config.ISL_CNT() * config.ISL_SIZE() == config.POP_SIZE());

      // set random pointer seed
      random_ptr = emp::NewPtr<emp::Random>(config.SEED());

      // initialize the world
      Initialize();
    }

    ~DiagWorld()
    {
      selection.Delete();
      random_ptr.Delete();
      diagnostic.Delete();
      pop_fit.Delete();
      pop_opti.Delete();
      pnt_fit.Delete();
      pnt_opti.Delete();
      pop_str.Delete();
      structure.Delete();
    }

    ///< functions called to setup the world

    // call all functions to initiallize the world
    void Initialize();

    // set OnUpdate function from World.h
    void SetOnUpdate();

    // set mutation operator from World.h
    void SetMutation();

    // set selction scheme
    void SetSelection();

    // set what to do when offspring is ready to go
    void SetOnOffspringReady();

    // set evaluation function
    void SetEvaluation();

    // set data tracking with data nodes
    void SetDataTracking();

    // set structure object if not well-mixed population
    void SetStructure();

    // populate the world with initial solutions
    void PopulateWorld();


    ///< principle steps during an evolutionary run

    // reset all data step
    void ResetData();

    // evaluation step
    void EvaluationStep();

    // population structure step
    void StructureStep();

    // selction step
    void SelectionStep();

    // reprodutive step
    void ReproductionStep();

    // record data step
    void RecordData();


    ///< selection scheme implementations

    void Truncation();

    void Tournament();

    void FitnessSharing();

    void EpsilonLexicase();

    void NonDominatedSorting();

    void NoveltySearch();


    ///< diagnostic function implementations

    void ExploitationRate();

    void OrderedExploitation();

    void MultiPathExploration();

    void ContradictoryObjectives();

    void MultiValleyCrossing();


    ///< data tracking

    void FindEverything();

    void CollectActivationCoverage();

    void CollectSatisfactoryCoverage();

    double IslandFitnessAverage();


    // size_t ActivationGeneOverlap(); might be useful for island overlaps?

    double MaxPopTrait();

    double MaxPopGene();

    ///< helper functions

    // create a matrix of popultion phenotype vectors
    fmatrix_t PopFitMat();


  private:
    // experiment configurations
    DiaConfig & config;
    enum Scheme {TRUNCATION=0,TOURNAMENT=1,LEXICASE=2};

    // vector holding population aggregate scores (by position id)
    phenotype_t fit_vec;
    // vector holding parent solutions selected by selection scheme
    ids_t parent_vec;
    // novelty minimum
    double pmin = 0.0;
    // generations since solution added to archive
    size_t archive_gens = 0;

    // evaluation lambda we set
    eval_t evaluation;
    // selection lambda we set
    sele_t select;


    // select.h var
    emp::Ptr<Selection> selection;
    // problem.h var
    emp::Ptr<Diagnostic> diagnostic;
    // structure.h var
    emp::Ptr<Structure> structure = nullptr;

    ///< data file & node related variables

    // file we are working with
    emp::DataFile data_file;
    // node to track population fitnesses
    nodef_t pop_fit;
    // node to track population opitmized count
    nodeo_t pop_opti;
    // node to track parent fitnesses
    nodef_t pnt_fit;
    // node to track parent optimized count
    nodeo_t pnt_opti;
    // node to track streak counts
    nodeo_t pop_str;
    // csv file to track best performing solutions
    std::ofstream elite_pheno_csv;
    std::ofstream elite_geno_csv;

    ///< data we are tracking during an evolutionary run

    // elite solution position
    size_t elite_pos;
    // optimal solution position
    size_t opti_pos;
    // population activation gene vector
    emp::vector<size_t> pop_acti_gene;
    // population satisfactory trait coverage
    emp::vector<size_t> pop_satis_trt;

    // multi-valley crossing data
    // valley peaks for each floored integer gene value
    const phenotype_t peaks = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  8.0,  9.0,
                                9.0, 11.0, 11.0, 11.0, 14.0, 14.0, 14.0, 14.0, 18.0, 18.0,
                                18.0, 18.0, 18.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 29.0,
                                29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 36.0, 36.0, 36.0, 36.0,
                                36.0, 36.0, 36.0, 36.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0,
                                44.0, 44.0, 44.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0,
                                53.0, 53.0, 53.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0,
                                63.0, 63.0, 63.0, 63.0, 74.0, 74.0, 74.0, 74.0, 74.0, 74.0,
                                74.0, 74.0, 74.0, 74.0, 74.0, 74.0, 86.0, 86.0, 86.0, 86.0,
                                86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 99.0};
    // where do the dips start?
    const double dips_start = 8.0;
    // where do dips end?
    const double dips_end = 99.9;

};

///< functions called to setup the world

void DiagWorld::Initialize()
{
  std::cout << "==========================================" << std::endl;
  std::cout << "BEGINNING INITIAL SETUP" << std::endl;
  std::cout << "==========================================" << std::endl;

  // reset the world upon start
  Reset();
  // set world to well mixed so we don't over populate
  SetPopStruct_Mixed(true);


  // stuff we need to initialize for the experiment
  SetEvaluation();
  SetMutation();
  SetOnUpdate();
  SetDataTracking();
  SetSelection();
  SetOnOffspringReady();
  SetStructure();
  PopulateWorld();

  std::cout << "==========================================" << std::endl;
  std::cout << "FINISHED INITIAL SETUP" << std::endl;
  std::cout << "==========================================" << std::endl;
}

void DiagWorld::SetOnUpdate()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting OnUpdate function..." << std::endl;

  // set up the evolutionary algorithm
  OnUpdate([this](size_t gen)
  {

    // step 0: reset all data collection variables
    ResetData();

    // step 1: evaluate all solutions on diagnostic
    EvaluationStep();

    // step 2: create structure if needed
    StructureStep();

    // step 3: select parent solutions for
    SelectionStep();

    // step 4: gather and record data
    RecordData();

    // step 5: reproduce and create new solutions
    ReproductionStep();
  });

  std::cout << "Finished setting the OnUpdate function! \n" << std::endl;
}

void DiagWorld::SetMutation()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting mutation function..." << std::endl;

  // set the mutation function
  SetMutFun([this](Org & org, emp::Random & random)
  {
    // number of mutations and solution genome
    size_t mcnt = 0;
    genome_t & genome = org.GetGenome();

    // quick checks
    emp_assert(genome.size() == config.DIMENSIONALITY());
    emp_assert(0.0 < config.TARGET());

    for(size_t i = 0; i < genome.size(); ++i)
    {
      // if we do a mutation at this objective
      if(random_ptr->P(config.MUTATE_PER()))
      {
        const double mut = random_ptr->GetRandNormal(config.MEAN(), config.STD());

        // rebound if
        if(config.TARGET() < genome[i] + mut)
        {
          genome[i] = config.TARGET() - (genome[i] + mut - config.TARGET());
        }
        // rebound if
        else if(genome[i] + mut < config.LOWER_BND())
        {
          genome[i] = std::abs(genome[i] + mut) + config.LOWER_BND();
        }
        // else add mutation
        else
        {
          genome[i] = genome[i] + mut;
        }
        ++mcnt;
      }
    }

    return mcnt;
  });

  std::cout << "Mutation function set!\n" << std::endl;
}

void DiagWorld::SetSelection()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting Selection function..." << std::endl;

  selection = emp::NewPtr<Selection>(random_ptr);
  std::cout << "Created selection" << std::endl;

  switch (config.SELECTION())
  {
    case static_cast<size_t>(Scheme::TRUNCATION):
      Truncation();
      break;

    case static_cast<size_t>(Scheme::TOURNAMENT):
      Tournament();
      break;

    case static_cast<size_t>(Scheme::LEXICASE):
      EpsilonLexicase();
      break;

    default:
      std::cout << "ERROR UNKNOWN SELECTION CALL" << std::endl;
      emp_assert(true);
      break;
  }

  std::cout << "Finished setting the Selection function! \n" << std::endl;
}

void DiagWorld::SetOnOffspringReady()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting OnOffspringReady function..." << std::endl;

  OnOffspringReady([this](Org & org, size_t parent_pos)
  {
    // quick checks
    emp_assert(fun_do_mutations); emp_assert(random_ptr);
    emp_assert(org.GetGenome().size() == config.DIMENSIONALITY());
    emp_assert(org.GetM() == config.DIMENSIONALITY());

    // do mutations on offspring
    size_t mcnt = fun_do_mutations(org, *random_ptr);

    // get parent
    Org & parent = *pop[parent_pos];

    // no mutations were applied to offspring
    if(mcnt == 0)
    {
      // quick checks
      emp_assert(parent.GetGenome().size() == config.DIMENSIONALITY());
      emp_assert(parent.GetM() == config.DIMENSIONALITY());

      // give everything to offspring from parent
      org.MeClone();
      org.Inherit(parent.GetPhenotype(), parent.GetSatisfactoryVec(), parent.GetCount(), parent.GetAggregate(), parent.GetStart(), parent.GetStreak(), parent.GetIslandID());
    }
    else
    {
      org.Reset();
      org.SetIslandID(parent.GetIslandID());
    }
  });

  std::cout << "Finished setting OnOffspringReady function!\n" << std::endl;
}

void DiagWorld::SetEvaluation()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting Evaluation function..." << std::endl;

  diagnostic = emp::NewPtr<Diagnostic>(config.TARGET(), config.LOWER_BND());
  std::cout << "Created diagnostic emp::Ptr" << std::endl;

  switch (config.DIAGNOSTIC())
  {
    case 0: // exploitation
      ExploitationRate();
      break;

    case 1: // structured exploitation
      OrderedExploitation();
      break;

    case 2: // contradictory objectives
      ContradictoryObjectives();
      break;

    case 3: // exploration
      MultiPathExploration();
      break;

    default: // error, unknown diganotic
      std::cout << "ERROR: UNKNOWN DIAGNOSTIC" << std::endl;
      emp_assert(false);
      break;
  }

  std::cout << "Evaluation function set!\n" <<std::endl;
}

void DiagWorld::SetDataTracking()
{
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Setting up data tracking..." << std::endl;

  // initialize all nodes
  std::cout << "Initializing nodes..." << std::endl;
  pop_fit.New();
  pop_opti.New();
  pnt_fit.New();
  pnt_opti.New();
  pop_str.New();
  std::cout << "Nodes initialized!" << std::endl;

  // update we are at
  data_file.AddFun<size_t>([this]()
  {
    return update;
  }, "gen", "Current generation at!");

  // track population aggregate phenotype stats: average, variance, min, max
  data_file.AddMean(*pop_fit, "pop_fit_avg", "Population average aggregate performance.");
  data_file.AddVariance(*pop_fit, "pop_fit_var", "Population variance aggregate performance.");
  data_file.AddMax(*pop_fit, "pop_fit_max", "Population maximum aggregate performance.");
  data_file.AddMin(*pop_fit, "pop_fit_min", "Population minimum aggregate performance.");

  // track population optimized objective count stats: average, variance, min, max
  data_file.AddMean(*pop_opti, "pop_opt_avg", "Population average objective optimization count.");
  data_file.AddVariance(*pop_opti, "pop_opt_var", "Population variance objective optimization count.");
  data_file.AddMax(*pop_opti, "pop_opt_max", "Population maximum objective optimization count.");
  data_file.AddMin(*pop_opti, "pop_opt_min", "Population minimum objective optimization count.");

  // track parent aggregate phenotype stats: average, variance, min, max
  data_file.AddMean(*pnt_fit, "pnt_fit_avg", "Parent average aggregate performance.");
  data_file.AddVariance(*pnt_fit, "pnt_fit_var", "Parent variance aggregate performance.");
  data_file.AddMax(*pnt_fit, "pnt_fit_max", "Parent maximum aggregate performance.");
  data_file.AddMin(*pnt_fit, "pnt_fit_min", "Parent minimum aggregate performance.");

  // track parent optimized objective count stats: average, variance, min, max
  data_file.AddMean(*pnt_opti, "pnt_opt_avg", "Parent average objective optimization count.");
  data_file.AddVariance(*pnt_opti, "pnt_opt_var", "Parent variance objective optimization count.");
  data_file.AddMax(*pnt_opti, "pnt_opt_max", "Parent maximum objective optimization count.");
  data_file.AddMin(*pnt_opti, "pnt_opt_min", "Parent minimum objective optimization count.");

  // track parent optimized objective count stats: average, variance, min, max
  data_file.AddMean(*pop_str, "pop_str_avg", "Population average streak count.");
  data_file.AddVariance(*pop_str, "pop_str_var", "Population variance streak count.");
  data_file.AddMax(*pop_str, "pop_str_max", "Population maximum streak count.");
  data_file.AddMin(*pop_str, "pop_str_min", "Population minimum streak count.");

  std::cout << "Added all data nodes to data file!" << std::endl;

  // population satisfactory traits count
  data_file.AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(pop_satis_trt.size() == config.DIMENSIONALITY());
    return std::count_if(pop_satis_trt.begin(), pop_satis_trt.end(), [](size_t i) { return 0 < i;});

  }, "pop_sat_cov", "Satisfactory trait coverage in the population!");

  // unique island satisfactory traits count
  data_file.AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(pop_satis_trt.size() == config.DIMENSIONALITY());
    return std::count_if(pop_satis_trt.begin(), pop_satis_trt.end(), [](size_t i) { return 1 == i;});

  }, "isl_sat_cov", "Unique satisfactory trait coverage in the population!");

  // population activation gene coverage
  data_file.AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(pop_acti_gene.size() == config.DIMENSIONALITY());
    return std::count_if(pop_acti_gene.begin(), pop_acti_gene.end(), [](size_t i) { return 0 < i;});

  }, "pop_act_cov", "Activation gene in the population!");

  // unique island activation gene coverage
  data_file.AddFun<size_t>([this]()
  {
    // quick checks
    emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(pop_acti_gene.size() == config.DIMENSIONALITY());
    return std::count_if(pop_acti_gene.begin(), pop_acti_gene.end(), [](size_t i) { return 1 == i;});

  }, "isl_act_cov", "Unique activation genes in the population !");

  // average fitness of best individual per island
  data_file.AddFun<double>([this]()
  {
    // quick checks
    emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(pop_acti_gene.size() == config.DIMENSIONALITY());
    return IslandFitnessAverage();

  }, "isl_fit_avg", "Average fitness from best individual per island !");

  // loss of diversity
  data_file.AddFun<double>([this]()
  {
    // quick checks
    emp_assert(parent_vec.size() == config.POP_SIZE());

    std::set<size_t> unique;
    for(auto & id : parent_vec) {unique.insert(id);}

    // ask Charles
    const double num = static_cast<double>(unique.size());
    const double dem = static_cast<double>(config.POP_SIZE());

    return num / dem;
  }, "los_div", "Loss in diversity generated by the selection scheme!");

  // selection pressure
  data_file.AddFun<double>([this]()
  {
    // quick checks
    emp_assert(pop_fit->GetCount() == config.POP_SIZE());
    emp_assert(pnt_fit->GetCount() == config.POP_SIZE());

    const double pop = pop_fit->GetMean();
    const double pnt = pnt_fit->GetMean();
    const double var = pop_fit->GetVariance();

    if(var == 0.0) {return 0.0;}

    return (pop - pnt) / var;
  }, "sel_pre", "Selection pressure applied by selection scheme!");

  // selection variance
  data_file.AddFun<double>([this]()
  {
    // quick checks
    emp_assert(pop_fit->GetCount() == config.POP_SIZE());
    emp_assert(pnt_fit->GetCount() == config.POP_SIZE());

    const double pop = pop_fit->GetVariance();
    const double pnt = pnt_fit->GetVariance();

    if(pnt == 0.0) {return 0.0;}

    return pop / pnt;
  }, "sel_var", "Selection pressure applied by selection scheme!");

  // max trait in the population
  data_file.AddFun<double>([this]()
  {
    return MaxPopTrait();
  }, "pop_max_trt", "Maximum trait value found in the population!");

  // max gene in the population
  data_file.AddFun<double>([this]()
  {
    return MaxPopGene();
  }, "pop_max_gene", "Maximum gene value found in the population!");

  data_file.PrintHeaderKeys();

  // create elite csv plus headers
  elite_pheno_csv.open(config.OUTPUT_DIR() + "elite-pheno.csv");
  elite_geno_csv.open(config.OUTPUT_DIR() + "elite-geno.csv");

  std::string header_pheno = "Gen";
  std::string header_geno = "Gen";
  for(size_t i = 0; i < config.DIMENSIONALITY(); ++i)
  {
    header_pheno += ",t";
    header_geno += ",g";

    header_pheno += std::to_string(i);
    header_geno += std::to_string(i);
  }

  elite_pheno_csv << header_pheno << "\n";
  elite_geno_csv << header_geno << "\n";

  std::cout << "Finished setting data tracking!\n" << std::endl;
}

void DiagWorld::SetStructure()
{
  structure = emp::NewPtr<Structure>(config.ISL_CNT(), config.ISL_SIZE(), random_ptr);

  if(config.ISL_CNT() == 1) {structure->GenerateStructure(ids_t(config.POP_SIZE(),0));}
}

void DiagWorld::PopulateWorld()
{
  // quick checks
  emp_assert(structure != nullptr);

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Populating world with initial solutions..." << std::endl;

  // initialize population with random vectors between LOWER_BND & UPPER_BND
  for(int i = 0; i < config.POP_SIZE(); ++i)
  {
    genome_t g = emp::RandomDoubleVector(*random_ptr, config.DIMENSIONALITY(), config.LOWER_BND(), config.UPPER_BND());
    Inject(g,1);
  }

  // quick checks
  emp_assert(0 < pop.size()); emp_assert(pop.size() == config.POP_SIZE());

  // add population structure for selection scheme

  // well-mixed
  if(config.ISL_CNT() == 1)
  {
    // same island for all
    for(size_t i = 0; i < pop.size(); ++i)
    {
      Org & org = *pop[i];
      org.SetIslandID(0);
    }
  }
  // if an island struture
  else
  {
    // will add solutions to their starting island
    for(size_t i = 0; i < config.POP_SIZE(); ++i)
    {
      Org & org = *pop[i];
      org.SetIslandID(i % config.ISL_CNT());
    }
  }

  std::cout << "Initialing world complete!" << std::endl;
}


///< principle steps during an evolutionary run

void DiagWorld::ResetData()
{
  // reset all data nodes
  pop_fit->Reset();
  pop_opti->Reset();
  pnt_fit->Reset();
  pnt_opti->Reset();
  pop_str->Reset();

  // reset all positon ids
  elite_pos = config.POP_SIZE();
  opti_pos = config.POP_SIZE();

  // reset all vectors/maps holding current gen data
  fit_vec.clear();
  parent_vec.clear();
  pop_acti_gene.clear();
  pop_satis_trt.clear();
}

void DiagWorld::EvaluationStep()
{
  // quick checks
  emp_assert(fit_vec.size() == 0); emp_assert(0 < pop.size());
  emp_assert(pop.size() == config.POP_SIZE());

  // iterate through the world and populate fitness vector
  fit_vec.resize(config.POP_SIZE());
  for(size_t i = 0; i < pop.size(); ++i)
  {
    Org & org = *pop[i];

    // no evaluation needed if offspring is a clone
    fit_vec[i] = (org.GetClone()) ? org.GetAggregate() : evaluation(org);
  }
}

void DiagWorld::StructureStep()
{
  // quick checks
  emp_assert(parent_vec.size() == 0); emp_assert(0 < pop.size());
  emp_assert(pop.size() == config.POP_SIZE()); emp_assert(structure != nullptr);

  // we only do something if the population is not well-mixed
  if(1 == config.ISL_CNT()) {return;}

  // do we need to execute a migration event?
  if(config.MIGRATION())
  {
    if((update % config.MIGRATION_INT() == 0) and (0 < update))
    {
      // migration time
      structure->Migration(config.MIGRATION_SZE());

      // go through each island in the structure
      for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
      {
        for(const size_t & id : structure->GetIsland(isl))
        {
          Org & org = *pop[id];
          org.SetIslandID(isl);
        }
      }
    }
  }

  // update the structure with the current set of org island ids
  ids_t isl_ids(config.POP_SIZE());
  for(size_t i = 0; i < config.POP_SIZE(); ++i)
  {
    const Org & org = *pop[i];
    isl_ids[i] = org.GetIslandID();
  }

  structure->GenerateStructure(isl_ids);
}

void DiagWorld::SelectionStep()
{
  // quick checks
  emp_assert(parent_vec.size() == 0); emp_assert(0 < pop.size());
  emp_assert(pop.size() == config.POP_SIZE());

  // store parents
  auto parents = select();
  emp_assert(parents.size() == config.POP_SIZE());

  parent_vec = parents;
}

void DiagWorld::RecordData()
{
  /// Add data to all nodes

  // get pop data
  emp_assert(pop.size() == config.POP_SIZE());
  for(size_t i = 0; i < pop.size(); ++i)
  {
    const Org & org = *pop[i];
    pop_fit->Add(org.GetAggregate());
    pop_opti->Add(org.GetCount());
    pop_str->Add(org.GetStreak());

    const Org & par = *pop[parent_vec[i]];
    pnt_fit->Add(par.GetAggregate());
    pnt_opti->Add(par.GetCount());
  }
  emp_assert(pop_fit->GetCount() == config.POP_SIZE());
  emp_assert(pop_opti->GetCount() == config.POP_SIZE());
  emp_assert(pop_str->GetCount() == config.POP_SIZE());

  // get parent data
  emp_assert(parent_vec.size() == config.POP_SIZE());
  emp_assert(pnt_fit->GetCount() == config.POP_SIZE());
  emp_assert(pnt_opti->GetCount() == config.POP_SIZE());

  /// get all position ids

  FindEverything();

  /// fill vectors & map
  emp_assert(fit_vec.size() == config.POP_SIZE()); // should be set already
  emp_assert(parent_vec.size() == config.POP_SIZE()); // should be set already

  /// update the file
  data_file.Update();

  // record elite solution traits
  Org & ele = *pop[elite_pos];

  std::string traits = std::to_string(update);
  const auto & p = ele.GetPhenotype();
  for(size_t i = 0; i < p.size(); ++i)
  {
    traits += ",";
    traits += std::to_string(p[i]);
  }
  elite_pheno_csv << traits << "\n";

  std::string genes = std::to_string(update);
  const auto & g = ele.GetGenome();
  for(size_t i = 0; i < g.size(); ++i)
  {
    genes += ",";
    genes += std::to_string(g[i]);
  }
  elite_geno_csv << genes << "\n";


  // output this so we know where we are in terms of generations and fitness
  Org & org = *pop[elite_pos];
  Org & opt = *pop[opti_pos];
  std::cout << "gen=" << GetUpdate() << ", max_fit=" << org.GetAggregate()  << ", max_opt=" << opt.GetCount() << std::endl;
}

void DiagWorld::ReproductionStep()
{
  // quick checks
  emp_assert(parent_vec.size() == config.POP_SIZE());
  emp_assert(pop.size() == config.POP_SIZE());

  // go through parent ids and do births
  for(auto & id : parent_vec)
  {
    DoBirth(GetGenomeAt(id), id);
  }
}


///< selection scheme implementations

void DiagWorld::Truncation()
{
  std::cout << "Setting selection scheme: Truncation" << std::endl;

  // set select using truncation selection
  select = [this]()
  {
    // quick checks
    emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
    emp_assert(0 < pop.size()); emp_assert(fit_vec.size() == config.POP_SIZE());
    emp_assert(structure != nullptr);

    if(config.ISL_CNT() == 1)
    {
      // group population by fitness
      const fitgp_t group = selection->FitnessGroup(fit_vec);
      return selection->Truncation(config.TRUNC_SIZE(), config.POP_SIZE(), group);
    }
    else
    {
      // do selection with island structure
      ids_t parents;

      for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
      {
        phenotype_t isl_fits;
        for(const size_t & id : structure->GetIsland(isl))
        {
          Org & org = *pop[id];
          isl_fits.push_back(org.GetAggregate());
        }
        emp_assert(isl_fits.size() == config.ISL_SIZE());

        // get parents for the island
        const fitgp_t group = selection->FitnessGroup(isl_fits);
        const ids_t isl_parents = selection->Truncation(config.TRUNC_SIZE(), config.ISL_SIZE(), group);
        emp_assert(isl_parents.size() == config.ISL_SIZE());

        for(const size_t & id : isl_parents) {parents.push_back(structure->GetIslandID(isl, id));}
      }
      emp_assert(parents.size() == config.POP_SIZE());

      return parents;
    }
  };

  std::cout << "Truncation selection scheme set!" << std::endl;
}

void DiagWorld::Tournament()
{
  std::cout << "Setting selection scheme: Tournament" << std::endl;

  select = [this]()
  {
    if(config.ISL_CNT() == 1)
    {
      // quick checks
      emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
      emp_assert(0 < pop.size()); emp_assert(fit_vec.size() == config.POP_SIZE());

      // will hold parent ids + get pop agg score values
      ids_t parent(pop.size());

      // get pop size amount of parents
      for(size_t i = 0; i < parent.size(); ++i)
      {
        parent[i] = selection->Tournament(config.TOUR_SIZE(), fit_vec);
      }

      return parent;
    }
    else
    {
      // quick checks
      emp_assert(structure != nullptr);

      // do selection with island structure
      ids_t parents;

      for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
      {
        // gather fitnesses on the island
        phenotype_t isl_fits;
        for(const size_t & id : structure->GetIsland(isl))
        {
          Org & org = *pop[id];
          isl_fits.push_back(org.GetAggregate());
        }
        emp_assert(isl_fits.size() == config.ISL_SIZE());

        // run tournament selection on each island
        for(size_t i = 0; i < config.ISL_SIZE(); ++i)
        {
          const size_t winner = selection->Tournament(config.TOUR_SIZE(), isl_fits);
          parents.push_back(structure->GetIslandID(isl, winner));
        }
        emp_assert(parents.size() == config.ISL_SIZE() * (isl + 1));
      }
      emp_assert(parents.size() == config.POP_SIZE());

      return parents;
    }
  };

  std::cout << "Tournament selection scheme set!" << std::endl;
}

void DiagWorld::EpsilonLexicase()
{
  std::cout << "Setting selection scheme: EpsilonLexicase" << std::endl;
  std::cout << "Epsilon: " << config.LEX_EPS() << std::endl;

  select = [this]()
  {
    if(config.ISL_CNT() == 1)
    {
      // quick checks
      emp_assert(selection); emp_assert(pop.size() == config.POP_SIZE());
      emp_assert(0 < pop.size());

      const fmatrix_t matrix = PopFitMat();

      // select parent ids
      ids_t parent(pop.size());

      for(size_t i = 0; i < parent.size(); ++i)
      {
        parent[i] = selection->EpsiLexicase(matrix, config.LEX_EPS(), config.DIMENSIONALITY());
      }

      return parent;
    }
    else
    {
      // quick checks
      emp_assert(structure != nullptr);

      // do selection with island structure
      ids_t parents;

      for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
      {
        // gather fitnesses on the island
        fmatrix_t isl_phenos;
        for(const size_t & id : structure->GetIsland(isl))
        {
          Org & org = *pop[id];
          isl_phenos.push_back(org.GetPhenotype());
        }
        emp_assert(isl_phenos.size() == config.ISL_SIZE());

        // run tournament selection on each island
        for(size_t i = 0; i < config.ISL_SIZE(); ++i)
        {
          const size_t winner = selection->EpsiLexicase(isl_phenos, config.LEX_EPS(), config.DIMENSIONALITY());
          parents.push_back(structure->GetIslandID(isl, winner));
        }
        emp_assert(parents.size() == config.ISL_SIZE() * (isl + 1));
      }
      emp_assert(parents.size() == config.POP_SIZE());

      return parents;
    }
  };

  std::cout << "Epsilon Lexicase selection scheme set!" << std::endl;
}

///< evaluation function implementations

void DiagWorld::ExploitationRate()
{
  std::cout << "Setting exploitation diagnostic..." << std::endl;

  evaluation = [this](Org & org)
  {
    // set phenotype & aggregate
    phenotype_t phenotype = diagnostic->ExploitationRate(org.GetGenome());

    //check if we are adding multi-valley crossing
    if(config.VALLEY_CROSSING()) {phenotype = diagnostic->MultiValleyCrossing(phenotype, peaks, dips_start, dips_end);}

    // set org evals
    org.SetPhenotype(phenotype);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->SatisfactoryVector(org.GetPhenotype(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Exploitation diagnotic set!" << std::endl;
}

void DiagWorld::OrderedExploitation()
{
  std::cout << "Setting structured exploitation diagnostic..." << std::endl;

  evaluation = [this](Org & org)
  {
    // set phenotype & aggregate
    phenotype_t phenotype = diagnostic->OrderedExploitation(org.GetGenome());

    //check if we are adding multi-valley crossing
    if(config.VALLEY_CROSSING()) {phenotype = diagnostic->MultiValleyCrossing(phenotype, peaks, dips_start, dips_end);}

    // set org evals
    org.SetPhenotype(phenotype);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->SatisfactoryVector(org.GetPhenotype(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Structured exploitation diagnotic set!" << std::endl;
}

void DiagWorld::MultiPathExploration()
{
  std::cout << "Setting exploration diagnostic..." << std::endl;

  evaluation = [this](Org & org)
  {
    // set phenotype & aggregate
    phenotype_t phenotype = diagnostic->MultiPathExploration(org.GetGenome());

    //check if we are adding multi-valley crossing
    if(config.VALLEY_CROSSING()) {phenotype = diagnostic->MultiValleyCrossing(phenotype, peaks, dips_start, dips_end);}

    // set org evals
    org.SetPhenotype(phenotype);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->SatisfactoryVector(org.GetPhenotype(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Exploration diagnotic set!" << std::endl;
}

void DiagWorld::ContradictoryObjectives()
{
  std::cout << "Setting contradictory objectives diagnostic..." << std::endl;

  evaluation = [this](Org & org)
  {
    // set phenotype & aggregate
    phenotype_t phenotype = diagnostic->ContradictoryObjectives(org.GetGenome());

    //check if we are adding multi-valley crossing
    if(config.VALLEY_CROSSING()) {phenotype = diagnostic->MultiValleyCrossing(phenotype, peaks, dips_start, dips_end);}

    // set org evals
    org.SetPhenotype(phenotype);
    org.AggregateScore();

    // set the starting position
    org.StartPosition();

    // set optimal vector and count
    optimal_t opti = diagnostic->SatisfactoryVector(org.GetPhenotype(), config.ACCURACY());
    org.SetOptimal(opti);
    org.CountOptimized();

    // set streak
    org.CalcStreak();

    return org.GetAggregate();
  };

  std::cout << "Weak ecology diagnotic set!" << std::endl;
}


///< data tracking

void DiagWorld::FindEverything()
{
  // quick checks
  emp_assert(fit_vec.size() == config.POP_SIZE()); emp_assert(pop_acti_gene.size() == 0);
  emp_assert(pop_satis_trt.size() == 0); emp_assert(elite_pos == config.POP_SIZE());

  // bools to make sure got everything
  bool elite_b = false,  opti_b = false;

  // collect coverages
  pop_acti_gene = ids_t(config.DIMENSIONALITY(), 0);
  pop_satis_trt = ids_t(config.DIMENSIONALITY(), 0);

  // loop and get data
  for(size_t i = 0; i < pop.size(); ++i)
  {
    const Org & org = *pop[i];

    // check if we need to do anything below
    if(elite_b and opti_b) {continue;}

    // find max fit solution
    if(org.GetAggregate() == pop_fit->GetMax() and !elite_b) {elite_b = true; elite_pos = i;}
    // find max satisfactory trait count solution
    if(org.GetCount() == pop_opti->GetMax() and !opti_b) {opti_b = true; opti_pos = i;}
  }

  // collect activation gene coverage
  CollectActivationCoverage();

  // collect satisfactory traits data
  CollectSatisfactoryCoverage();
}

void DiagWorld::CollectActivationCoverage()
{
  // collect activation gene data
  for(size_t obj = 0; obj < config.DIMENSIONALITY(); ++obj)
  {
    for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
      {
        for(const size_t & id : structure->GetIsland(isl))
        {
          const Org & org = *pop[id];
          if(org.GetStart() == obj)
          {
            pop_acti_gene[org.GetStart()]++;
            break;
          }
        }
      }
  }
}

void DiagWorld::CollectSatisfactoryCoverage()
{
  // collect satisfactory traits data
  for(size_t obj = 0; obj < config.DIMENSIONALITY(); ++obj)
  {
    for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
      {
        for(const size_t & id : structure->GetIsland(isl))
        {
          const Org & org = *pop[id];
          emp_assert(org.GetSatisfactoryVec().size() == config.DIMENSIONALITY());

          if(org.OptimizedAt(obj))
          {
            pop_satis_trt[obj]++;
            break;
          }
        }
      }
  }
}

double DiagWorld::IslandFitnessAverage()
{
  // quick checks
  emp_assert(fit_vec.size() == config.POP_SIZE()); emp_assert(structure != nullptr);

  // if well-mixed
  if(config.ISL_CNT() == 1)
  {
    return pop_fit->GetMax();
  }

  // best fitness per island
  emp::vector<double> isl_best_fit(config.ISL_CNT(), 0.0);

  // find best fitness per island
  for(size_t isl = 0; isl < config.ISL_CNT(); ++isl)
  {
    for(const size_t & id : structure->GetIsland(isl))
    {
      if(isl_best_fit[isl] < fit_vec[id]) {isl_best_fit[isl] = fit_vec[id];}
    }
  }

  return std::accumulate(isl_best_fit.begin(), isl_best_fit.end(), 0.0) / static_cast<double>(config.ISL_CNT());
}


double DiagWorld::MaxPopTrait()
{
  // iterate pop to check is a solution has the objective optimized
  double max = -1000000.0;
  for(size_t p = 0; p < pop.size(); ++p)
  {
    Org & org = *pop[p];

    if(max < org.GetMaxTrait())
    {
      max = org.GetMaxTrait();
    }
  }

  return max;
}

double DiagWorld::MaxPopGene()
{
  // iterate pop to check is a solution has the objective optimized
  double max = 0.0;
  for(size_t p = 0; p < pop.size(); ++p)
  {
    Org & org = *pop[p];

    for(auto & g : org.GetGenome())
    {
      if(max < g) {max = g;}
    }
  }

  return max;
}

///< helper functions

DiagWorld::fmatrix_t DiagWorld::PopFitMat()
{
  // quick checks
  emp_assert(pop.size() == config.POP_SIZE());

  // create matrix of population phenotype vectors
  fmatrix_t matrix(pop.size());

  for(size_t i = 0; i < pop.size(); ++i)
  {
    Org & org = *pop[i];
    emp_assert(org.GetPhenotype().size() == config.DIMENSIONALITY());

    // charles ask if this is the actual org phenotype vector or a deep copy made
    matrix[i] = org.GetPhenotype();
  }

  return matrix;
}

#endif