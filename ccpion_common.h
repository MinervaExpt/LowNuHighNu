#ifndef playlist_methods_h
#define playlist_methods_h

#include <algorithm>  // transform
#include <string>

#include "TString.h"  // Form

//RDF// This is eventually going to point to p4 MAD tuples and can be shared b/w
//RDF// me and Ben/Ever
std::string GetPlaylistFile(std::string plist, bool is_mc,
                            bool use_xrootd = true) {
  // const std::string processing_date = "20200713"; // new short tracking branches
  const std::string processing_date = "20211115";  // new recoil energy branches
  //  const std::string processing_date = "test"; // For test with small MAD tuplas
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir =
      is_mc ? "/minerva/data/users/granados/MAD_ana_plists/"
            : "/minerva/data/users/granados/MAD_ana_plists/";  // correct merging method
  topdir += processing_date;
  std::string playlist_file =
      use_xrootd ? Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(),
                        is_mc_str.c_str(), plist.c_str())
                 : Form("%s/%s_%s_plist.txt", topdir.c_str(), is_mc_str.c_str(),
                        plist.c_str());
  return playlist_file;
}

std::string GetTestPlaylist(bool is_mc) {
  return is_mc ? "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/ME1A_mc_p3_beta.txt" :
                 "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/ME1A_data_p3_beta.txt";
}

#endif  // playlist_methods_h
