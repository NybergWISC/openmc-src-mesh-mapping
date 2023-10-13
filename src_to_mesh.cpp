

#ifdef SRC_MAP_OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "openmc/capi.h"
#include "openmc/mesh.h"
#include "openmc/progress_bar.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/timer.h"

int main(int argc, char** argv) {

  openmc::settings::verbosity = 0;

  openmc_init(0, argv, nullptr);

  int total_sites = 1e6;
  int max_sites_incr = 10000;
  int sites_sampled {};

  #ifdef SRC_MAP_OPENMP
  std::cout << "Number of threads used: " << omp_get_max_threads() << std::endl;
  #else
  std::cout << "Threading disabled in this compilation" << std::endl;
  #endif

  // accessing mesh
  const auto& mesh = openmc::model::meshes[openmc::model::mesh_map[1]];

  std::cout << "Mesh ID: " << mesh->id() << std::endl;
  std::cout << "Mesh bins: " << mesh->n_bins() << std::endl;

  std::vector<size_t> bin_counts(mesh->n_bins(), 0);

  std::vector<uint64_t> seeds;
  #ifdef SRC_MAP_OPENMP
    seeds.resize(omp_get_max_threads());
  #else
    seeds.resize(1);
  #endif

  for (int i = 0; i < seeds.size(); i++) seeds[i] = i;

  std::cout << "Sampling source sites and mapping to the mesh...\n";

  openmc::Timer timer;

  timer.start();

  ProgressBar prog_bar{};

  while (sites_sampled < total_sites) {

    int n_sites = std::min(total_sites - sites_sampled, max_sites_incr);

    // sample source sites
    std::vector<openmc::SourceSite> sites(n_sites);

#pragma omp parallel
{
    #ifdef SRC_MAP_OPENMP
    uint64_t* seed = &seeds[omp_get_thread_num()];
    #else
    uint64_t* seed = &seeds[0];
    #endif

    #pragma omp for
    for (int i = 0; i < sites.size(); i++) {
      sites[i] = openmc::sample_external_source(seed);
    }
}
    sites_sampled += n_sites;

    #pragma omp parallel for shared(bin_counts)
    for (int i = 0; i < sites.size(); i++) {
      bin_counts[mesh->get_bin(sites[i].r)] += 1;
    }

    prog_bar.set_value(100. * (double)sites_sampled / (double)total_sites);
  }

  timer.stop();

  double mapping_time = timer.elapsed();
  timer.reset();


  std::vector<double> rel_strengths(bin_counts.begin(), bin_counts.end());
  for (auto& rs : rel_strengths) { rs /= sites_sampled; }

  // Determine total source strength
  double total_strength = 0.0;
  for (auto& s : openmc::model::external_sources)
    total_strength += s->strength();

  std::cout << "Total external source strength: " << total_strength << std::endl;

  // TODO: extend to all mesh types
  openmc::UnstructuredMesh* umesh_ptr = dynamic_cast<openmc::UnstructuredMesh*>(mesh.get());
  if (!umesh_ptr) { std::cerr << "Non-unstructured mesh used" << std::endl; return 1; }

  std::cout << "Writing output file..." << std::endl;
  timer.start();

  // generate a text output as CSV
  // Each row contains: OpenMC bin #, centroid (xyz), volume, rel. src strength
  std::ofstream output("mesh_src_strengths.csv");
  output << "bin, cx, cy, cz, volume, rel. src" << std::endl;
  for (int i = 0; i < mesh->n_bins(); i++) {
    std::stringstream line {};
    line << i << ", ";
    auto centroid = umesh_ptr->centroid(i);
    line << centroid.x << ", " << centroid.y << ", " << centroid.z << ", ";
    line << mesh->volume(i) << ", ";
    line << rel_strengths[i] << std::endl;

    output << line.str();

  }
  output.close();

  timer.stop();
  double output_time = timer.elapsed();

  // report time
  std::cout << "--------------" << std::endl;
  std::cout << "Time Summary: " << std::endl;
  std::cout << "--------------" << std::endl;

  std::cout << "Source mapping (" << total_sites << " sites): " << mapping_time << "s" << std::endl;
  std::cout << "File output: " << output_time << " s" << std::endl;

  openmc_finalize();

  return 0;
}
