#ifndef runEffPurTable_C
#define runEffPurTable_C

#include <chrono>

#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"  // typedefs EventCount, t12 globals for timing
#include "includes/Cuts.h"
#include "includes/EventSelectionTable.h"
#include "includes/MacroUtil.h"
#include "playlist_methods.h"  // GetPlaylistFile

//==============================================================================
// Loop and fill
//==============================================================================
std::tuple<EventCount, EventCount> FillCounters(
    const CCPi::MacroUtil& util, CVUniverse* universe, const EDataMCTruth& type,
    const EventCount& s, const EventCount& b = EventCount()) {
  EventCount signal = s;
  EventCount bg = b;
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    std::tie(signal, bg) =
        ccpi_event::FillCounters(event, signal, bg);  // Does a lot of work
  }                                                   // events
  std::cout << "*** Done ***\n\n";
  return {signal, bg};
}

//==============================================================================
// Main
//==============================================================================
void runEffPurTable(int signal_definition_int = 0, const char* plist = "ALL") {
  auto start = std::chrono::steady_clock::now();
  // INIT MACRO UTILITY OBJECT
  const std::string macro("runEffPurTable");
  bool do_data = true, do_mc = true, do_truth = true;
  bool do_systematics = false, is_grid = false;
  bool use_xrootd = false;
  // std::string data_file_list = GetPlaylistFile(plist, false, use_xrootd);
  // std::string mc_file_list = GetPlaylistFile(plist, true, use_xrootd);
  bool is_mc = false;
  std::string data_file_list = GetTestPlaylist(is_mc);
  is_mc = true;
  std::string mc_file_list = GetTestPlaylist(is_mc);
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // EFFICIENCY/PURITY COUNTERS
  // typdef EventCount map<ECut, double>
  EventCount n_remaining_sig, n_remaining_bg, n_remaining_data;

  std::tie(n_remaining_data, std::ignore) =
      FillCounters(util, util.m_data_universe, kData, n_remaining_data);

  std::tie(n_remaining_sig, n_remaining_bg) =
      FillCounters(util, util.m_error_bands.at("cv").at(0), kMC,
                   n_remaining_sig, n_remaining_bg);

  std::tie(n_remaining_sig, n_remaining_bg) =
      FillCounters(util, util.m_error_bands_truth.at("cv").at(0), kTruth,
                   n_remaining_sig, n_remaining_bg);

  PrintEffPurTable(n_remaining_sig, n_remaining_bg, n_remaining_data,
                   util.m_data_pot, util.m_mc_pot);
  auto stop = std::chrono::steady_clock::now();
  std::cout
      << "running time: "
      << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
      << "sec\n";
  std::cout << "cluster creation time: " << t01.count() << "sec\n";
  std::cout << "cluster loop 1: " << t12.count() << "sec\n";
  std::cout << "cluster loop 2: " << t23.count() << "sec\n";
}

#endif
