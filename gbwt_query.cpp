#include <iostream>
#include <fstream>

#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>

#include <gbwtgraph/path_cover.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

const std::string tool_name = "GFA to GBWTGraph";

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  // Parse command line options.
  int c = 0;
  while((c = getopt(argc, argv, "")) != -1)
  {
    switch(c)
    {
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Check command line options.
  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  std::string base_name = argv[optind]; optind++;

  // Initial output.
  Version::print(std::cout, tool_name);
  gbwt::printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = gbwt::readTimer();

  std::cout << "Parsing GFA and building GBWT..." << std::endl;
  auto results = gfa_to_gbwt(base_name + GFA_EXTENSION);
  auto index = *(results.first);

  if(results.first.get() == nullptr || results.second.get() == nullptr)
  {
    std::cerr << "gfa2gbwt: Construction failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Serializing GBWT..." << std::endl;
  {
    std::string gbwt_name = base_name + gbwt::GBWT::EXTENSION;
    std::ofstream out(gbwt_name, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << gbwt_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    results.first->serialize(out);
    out.close();
  }

  std::cout << "Building GBWTGraph..." << std::endl;
  GBWTGraph graph(*(results.first), *(results.second));

  std::cout << "Serializing GBWTGraph..." << std::endl;
  {
    std::string graph_name = base_name + GBWTGraph::EXTENSION;
    std::ofstream out(graph_name, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << graph_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    graph.serialize(out);
    out.close();
  }
  std::cout << std::endl;

  gbwt::printStatistics(*(results.first), base_name);

  double seconds = gbwt::readTimer();

  std::cout << "GBWTGraph built in " << seconds - start << " seconds" << std::endl;
  std::cout << "Memory usage " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;


  size_t paths_per_component = 100000;
 
  std::vector<gbwt::node_type> frequent_path
  {
    gbwt::Node::encode(29, false),
    gbwt::Node::encode(2, false),
    gbwt::Node::encode(28, false)
  };
  std::vector<gbwt::node_type> rare_path
  {
    gbwt::Node::encode(29, false),
    gbwt::Node::encode(2, false)
  };

  //gbwt::GBWT cover = local_haplotypes(graph, index, paths_per_component, 3);
  //gbwt::GBWT cover = path_cover_gbwt(graph, paths_per_component, 3);
  gbwt::SearchState frequent_state = index.find(frequent_path.begin(), frequent_path.end()); 
  std::cout << "3 node_count: " << frequent_state.size() << std::endl; 
  double seconds2 = gbwt::readTimer();

  std::cout << "GBWTGraph query in " << seconds2 - seconds << " seconds" << std::endl;

  //gbwt::GBWT cover2 = local_haplotypes(graph, index, paths_per_component, 2);
  gbwt::SearchState rare_state = index.find(rare_path.begin(), rare_path.end()); 
  std::cout << "2 node_count: " << rare_state.size() << std::endl; 
  double seconds3 = gbwt::readTimer();

  std::cout << "GBWTGraph query in " << seconds3 - seconds2 << " seconds" << std::endl;
  
  gbwt::GBWT cover2 = path_cover_gbwt(graph, paths_per_component, 2);
  std::cout << "path_size: " << cover2.metadata.paths() << std::endl; 
  std::cout << "sample_size: " << cover2.metadata.samples() << std::endl; 
  std::cout << "haplo_size: " << cover2.metadata.haplotypes() << std::endl; 
  
  

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: gfa2gbwt base_name" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
