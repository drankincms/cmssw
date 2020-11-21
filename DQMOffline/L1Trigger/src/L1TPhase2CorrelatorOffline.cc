/**
 *  @file     L1TPhase2CorrelatorOffline.cc
 *  @authors  Dylan Rankin (MIT)
 *  @date     20/10/2020
 *  @version  0.1
 *
 */

#include "DQMOffline/L1Trigger/interface/L1TPhase2CorrelatorOffline.h"

#include "TLorentzVector.h"
#include "TGraph.h"

#include <iomanip>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace reco;
using namespace trigger;
using namespace edm;
using namespace std;

const std::map<std::string, unsigned int> L1TPhase2CorrelatorOffline::PlotConfigNames = {
    {"resVsPt", PlotConfig::resVsPt},
    {"resVsEta", PlotConfig::resVsEta},
    {"ptDist", PlotConfig::ptDist},
    {"etaDist", PlotConfig::etaDist}};

//
// -------------------------------------- Constructor --------------------------------------------
//
L1TPhase2CorrelatorOffline::L1TPhase2CorrelatorOffline(const edm::ParameterSet& ps)
    : phase2PFToken_(
          consumes<std::vector<l1t::PFCandidate>>(ps.getUntrackedParameter<edm::InputTag>("phase2L1PfInputTag"))),
      phase2PuppiToken_(
          consumes<std::vector<l1t::PFCandidate>>(ps.getUntrackedParameter<edm::InputTag>("phase2L1PuppiInputTag"))),
      genJetToken_(consumes<std::vector<reco::GenJet>>(ps.getUntrackedParameter<edm::InputTag>("genJetsInputTag"))),
      genParticleToken_(
          consumes<std::vector<reco::GenParticle>>(ps.getUntrackedParameter<edm::InputTag>("genParticlesInputTag"))),
      objs_(ps.getParameter<edm::ParameterSet>("objects")),
      isParticleGun_(ps.getParameter<bool>("isParticleGun")),
      histFolder_(ps.getParameter<std::string>("histFolder")),
      respresolFolder_(histFolder_ + "/respresol_raw"),
      histDefinitions_(dqmoffline::l1t::readHistDefinitions(ps.getParameterSet("histDefinitions"), PlotConfigNames)),
      h_L1PF_pt_(),
      h_L1PF_eta_(),
      h_L1Puppi_pt_(),
      h_L1Puppi_eta_(),
      h_L1PF_pt_mu_(),
      h_L1PF_eta_mu_(),
      h_L1Puppi_pt_mu_(),
      h_L1Puppi_eta_mu_(),
      h_L1PF_pt_el_(),
      h_L1PF_eta_el_(),
      h_L1Puppi_pt_el_(),
      h_L1Puppi_eta_el_(),
      h_L1PF_pt_pho_(),
      h_L1PF_eta_pho_(),
      h_L1Puppi_pt_pho_(),
      h_L1Puppi_eta_pho_(),
      h_L1PF_pt_ch_(),
      h_L1PF_eta_ch_(),
      h_L1Puppi_pt_ch_(),
      h_L1Puppi_eta_ch_(),
      h_L1PF_pt_nh_(),
      h_L1PF_eta_nh_(),
      h_L1Puppi_pt_nh_(),
      h_L1Puppi_eta_nh_(),
      h_L1PF_electron_ptratio_0p2_vs_pt_barrel_(),
      h_L1PF_electron_ptratio_0p2_vs_pt_endcap_(),
      h_L1PF_electron_ptratio_0p2_vs_pt_ecnotk_(),
      h_L1PF_electron_ptratio_0p2_vs_pt_hf_(),
      h_L1PF_electron_ptratio_0p2_vs_eta_(),
      h_L1Puppi_electron_ptratio_0p2_vs_pt_barrel_(),
      h_L1Puppi_electron_ptratio_0p2_vs_pt_endcap_(),
      h_L1Puppi_electron_ptratio_0p2_vs_pt_ecnotk_(),
      h_L1Puppi_electron_ptratio_0p2_vs_pt_hf_(),
      h_L1Puppi_electron_ptratio_0p2_vs_eta_(),
      h_L1PF_electron_ptratio_best_vs_pt_barrel_(),
      h_L1PF_electron_ptratio_best_vs_pt_endcap_(),
      h_L1PF_electron_ptratio_best_vs_pt_ecnotk_(),
      h_L1PF_electron_ptratio_best_vs_pt_hf_(),
      h_L1PF_electron_ptratio_best_vs_eta_(),
      h_L1Puppi_electron_ptratio_best_vs_pt_barrel_(),
      h_L1Puppi_electron_ptratio_best_vs_pt_endcap_(),
      h_L1Puppi_electron_ptratio_best_vs_pt_ecnotk_(),
      h_L1Puppi_electron_ptratio_best_vs_pt_hf_(),
      h_L1Puppi_electron_ptratio_best_vs_eta_(),
      h_L1PF_pizero_ptratio_0p2_vs_pt_barrel_(),
      h_L1PF_pizero_ptratio_0p2_vs_pt_endcap_(),
      h_L1PF_pizero_ptratio_0p2_vs_pt_ecnotk_(),
      h_L1PF_pizero_ptratio_0p2_vs_pt_hf_(),
      h_L1PF_pizero_ptratio_0p2_vs_eta_(),
      h_L1Puppi_pizero_ptratio_0p2_vs_pt_barrel_(),
      h_L1Puppi_pizero_ptratio_0p2_vs_pt_endcap_(),
      h_L1Puppi_pizero_ptratio_0p2_vs_pt_ecnotk_(),
      h_L1Puppi_pizero_ptratio_0p2_vs_pt_hf_(),
      h_L1Puppi_pizero_ptratio_0p2_vs_eta_(),
      h_L1PF_pizero_ptratio_best_vs_pt_barrel_(),
      h_L1PF_pizero_ptratio_best_vs_pt_endcap_(),
      h_L1PF_pizero_ptratio_best_vs_pt_ecnotk_(),
      h_L1PF_pizero_ptratio_best_vs_pt_hf_(),
      h_L1PF_pizero_ptratio_best_vs_eta_(),
      h_L1Puppi_pizero_ptratio_best_vs_pt_barrel_(),
      h_L1Puppi_pizero_ptratio_best_vs_pt_endcap_(),
      h_L1Puppi_pizero_ptratio_best_vs_pt_ecnotk_(),
      h_L1Puppi_pizero_ptratio_best_vs_pt_hf_(),
      h_L1Puppi_pizero_ptratio_best_vs_eta_(),
      h_L1PF_pion_ptratio_0p2_vs_pt_barrel_(),
      h_L1PF_pion_ptratio_0p2_vs_pt_endcap_(),
      h_L1PF_pion_ptratio_0p2_vs_pt_ecnotk_(),
      h_L1PF_pion_ptratio_0p2_vs_pt_hf_(),
      h_L1PF_pion_ptratio_0p2_vs_eta_(),
      h_L1Puppi_pion_ptratio_0p2_vs_pt_barrel_(),
      h_L1Puppi_pion_ptratio_0p2_vs_pt_endcap_(),
      h_L1Puppi_pion_ptratio_0p2_vs_pt_ecnotk_(),
      h_L1Puppi_pion_ptratio_0p2_vs_pt_hf_(),
      h_L1Puppi_pion_ptratio_0p2_vs_eta_(),
      h_L1PF_pion_ptratio_best_vs_pt_barrel_(),
      h_L1PF_pion_ptratio_best_vs_pt_endcap_(),
      h_L1PF_pion_ptratio_best_vs_pt_ecnotk_(),
      h_L1PF_pion_ptratio_best_vs_pt_hf_(),
      h_L1PF_pion_ptratio_best_vs_eta_(),
      h_L1Puppi_pion_ptratio_best_vs_pt_barrel_(),
      h_L1Puppi_pion_ptratio_best_vs_pt_endcap_(),
      h_L1Puppi_pion_ptratio_best_vs_pt_ecnotk_(),
      h_L1Puppi_pion_ptratio_best_vs_pt_hf_(),
      h_L1Puppi_pion_ptratio_best_vs_eta_(),
      h_L1PF_jet_ptratio_vs_pt_barrel_(),
      h_L1PF_jet_ptratio_vs_pt_endcap_(),
      h_L1PF_jet_ptratio_vs_pt_ecnotk_(),
      h_L1PF_jet_ptratio_vs_pt_hf_(),
      h_L1PF_jet_ptratio_vs_eta_(),
      h_L1Puppi_jet_ptratio_vs_pt_barrel_(),
      h_L1Puppi_jet_ptratio_vs_pt_endcap_(),
      h_L1Puppi_jet_ptratio_vs_pt_ecnotk_(),
      h_L1Puppi_jet_ptratio_vs_pt_hf_(),
      h_L1Puppi_jet_ptratio_vs_eta_() {
  edm::LogInfo("L1TPhase2CorrelatorOffline") << "Constructor "
                                             << "L1TPhase2CorrelatorOffline::L1TPhase2CorrelatorOffline " << std::endl;

  auto reconames = objs_.getParameterNamesForType<std::vector<edm::InputTag>>();
  for (const std::string& name : reconames) {
    reco_.emplace_back(L1TPhase2CorrelatorOffline::MultiCollection(objs_, name, consumesCollector()), RecoVars());
  }
}

//
// -- Destructor
//
L1TPhase2CorrelatorOffline::~L1TPhase2CorrelatorOffline() {
  edm::LogInfo("L1TPhase2CorrelatorOffline")
      << "Destructor L1TPhase2CorrelatorOffline::~L1TPhase2CorrelatorOffline " << std::endl;
}

//
// -------------------------------------- beginRun --------------------------------------------
//
void L1TPhase2CorrelatorOffline::dqmBeginRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
  edm::LogInfo("L1TPhase2CorrelatorOffline") << "L1TPhase2CorrelatorOffline::beginRun" << std::endl;

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  bZ_ = magneticField->inTesla(GlobalPoint(0, 0, 0)).z();
}

//
// -------------------------------------- bookHistos --------------------------------------------
//
void L1TPhase2CorrelatorOffline::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const&, edm::EventSetup const&) {
  edm::LogInfo("L1TPhase2CorrelatorOffline") << "L1TPhase2CorrelatorOffline::bookHistograms" << std::endl;

  // book at beginRun
  bookPhase2CorrelatorHistos(ibooker);
}

//
// -------------------------------------- Analyze --------------------------------------------
//
void L1TPhase2CorrelatorOffline::analyze(edm::Event const& e, edm::EventSetup const& eSetup) {
  edm::Handle<std::vector<reco::GenJet>> genjets;
  edm::Handle<std::vector<reco::GenParticle>> genparticles;
  e.getByToken(genJetToken_, genjets);
  e.getByToken(genParticleToken_, genparticles);

  std::vector<const reco::GenParticle*> prompts, taus;
  for (const reco::GenParticle& gen : *genparticles) {
    if (isParticleGun_) {
      if (gen.statusFlags().isPrompt() == 1)
        prompts.push_back(&gen);
      continue;
    }
    if ((gen.isPromptFinalState() || gen.isDirectPromptTauDecayProductFinalState()) &&
        (std::abs(gen.pdgId()) == 11 || std::abs(gen.pdgId()) == 13) && gen.pt() > 5) {
      prompts.push_back(&gen);
    } else if (gen.isPromptFinalState() && std::abs(gen.pdgId()) == 22 && gen.pt() > 10) {
      prompts.push_back(&gen);
    } else if (abs(gen.pdgId()) == 15 && gen.isPromptDecayed()) {
      taus.push_back(&gen);
    }
  }

  for (auto& recopair : reco_) {
    recopair.first.get(e);
  }

  edm::Handle<std::vector<l1t::PFCandidate>> l1pfs;
  edm::Handle<std::vector<l1t::PFCandidate>> l1pups;
  e.getByToken(phase2PFToken_, l1pfs);
  e.getByToken(phase2PuppiToken_, l1pups);

  if (!l1pfs.isValid() or !l1pups.isValid()) {
    edm::LogWarning("L1TPhase2CorrelatorOffline") << "invalid collection: vector<l1t::PFCandidate> " << std::endl;
    return;
  }

  for (const reco::GenJet& j1 : *genjets) {
    bool ok = true;
    const reco::Candidate* match = nullptr;
    for (const reco::GenParticle* gp : prompts) {
      if (::deltaR2(*gp, j1) < 0.16f) {
        if (match != nullptr) {
          ok = false;
          break;
        } else {
          match = gp;
        }
      }
    }
    if (!ok)
      continue;
    if (!match) {
      // look for a tau
      for (const reco::GenParticle* gp : taus) {
        if (::deltaR2(*gp, j1) < 0.16f) {
          if (match != nullptr) {
            ok = false;
            break;
          } else {
            match = gp;
          }
        }
      }
      if (!ok)
        continue;
      if (match != nullptr && match->numberOfDaughters() == 2 &&
          std::abs(match->daughter(0)->pdgId()) + std::abs(match->daughter(1)->pdgId()) == 211 + 16) {
        // one-prong tau, consider it a pion
        match = (std::abs(match->daughter(0)->pdgId()) == 211 ? match->daughter(0) : match->daughter(1));
      }
    }
    if (match != nullptr) {
      if (std::abs(match->pdgId()) == 15) {
        reco::Particle::LorentzVector pvis;
        for (unsigned int i = 0, n = match->numberOfDaughters(); i < n; ++i) {
          const reco::Candidate* dau = match->daughter(i);
          if (std::abs(dau->pdgId()) == 12 || std::abs(dau->pdgId()) == 14 || std::abs(dau->pdgId()) == 16) {
            continue;
          }
          pvis += dau->p4();
        }
        mc_.pt = pvis.Pt();
        mc_.eta = pvis.Eta();
        mc_.phi = pvis.Phi();
      } else {
        mc_.fillP4(*match);
        mc_.fillPropagated(*match, bZ_);
      }
      mc_.id = std::abs(match->pdgId());
      mc_.iso04 = j1.pt() / mc_.pt - 1;
      mc_.iso02 = 0;
      for (const auto& dptr : j1.daughterPtrVector()) {
        if (::deltaR2(*dptr, *match) < 0.04f) {
          mc_.iso02 += dptr->pt();
        }
      }
      mc_.iso02 = mc_.iso02 / mc_.pt - 1;
    } else {
      if (j1.pt() < 20)
        continue;
      mc_.fillP4(j1);
      mc_.id = 0;
      mc_.iso02 = 0;
      mc_.iso04 = 0;
    }
    mc_.iso08 = mc_.iso04;
    for (const reco::GenJet& j2 : *genjets) {
      if (&j2 == &j1)
        continue;
      if (::deltaR2(j1, j2) < 0.64f)
        mc_.iso08 += j2.pt() / mc_.pt;
    }
    for (auto& recopair : reco_) {
      recopair.second.fill(recopair.first.objects(),
                           recopair.first.prop() ? mc_.caloeta : mc_.eta,
                           recopair.first.prop() ? mc_.calophi : mc_.phi);

      if (abs(mc_.id) == 11 && fabs(mc_.iso04) < 0.05) {
        if (recopair.first.name() == "L1PF") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1PF_electron_ptratio_0p2_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_electron_ptratio_best_vs_pt_barrel_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1PF_electron_ptratio_0p2_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_electron_ptratio_best_vs_pt_endcap_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1PF_electron_ptratio_0p2_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_electron_ptratio_best_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1PF_electron_ptratio_0p2_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_electron_ptratio_best_vs_pt_hf_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          }
          h_L1PF_electron_ptratio_0p2_vs_eta_->Fill(mc_.eta, recopair.second.pt02 / mc_.pt);
          h_L1PF_electron_ptratio_best_vs_eta_->Fill(mc_.eta, recopair.second.ptbest / mc_.pt);
        }
        if (recopair.first.name() == "L1Puppi") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1Puppi_electron_ptratio_0p2_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_electron_ptratio_best_vs_pt_barrel_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1Puppi_electron_ptratio_0p2_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_electron_ptratio_best_vs_pt_endcap_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1Puppi_electron_ptratio_0p2_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_electron_ptratio_best_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1Puppi_electron_ptratio_0p2_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_electron_ptratio_best_vs_pt_hf_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          }
          h_L1Puppi_electron_ptratio_0p2_vs_eta_->Fill(mc_.eta, recopair.second.pt02 / mc_.pt);
          h_L1Puppi_electron_ptratio_best_vs_eta_->Fill(mc_.eta, recopair.second.ptbest / mc_.pt);
        }
      }
      if (abs(mc_.id) == 111 && fabs(mc_.iso04) < 0.05) {
        if (recopair.first.name() == "L1PF") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1PF_pizero_ptratio_0p2_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pizero_ptratio_best_vs_pt_barrel_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1PF_pizero_ptratio_0p2_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pizero_ptratio_best_vs_pt_endcap_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1PF_pizero_ptratio_0p2_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pizero_ptratio_best_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1PF_pizero_ptratio_0p2_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pizero_ptratio_best_vs_pt_hf_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          }
          h_L1PF_pizero_ptratio_0p2_vs_eta_->Fill(mc_.eta, recopair.second.pt02 / mc_.pt);
          h_L1PF_pizero_ptratio_best_vs_eta_->Fill(mc_.eta, recopair.second.ptbest / mc_.pt);
        }
        if (recopair.first.name() == "L1Puppi") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1Puppi_pizero_ptratio_0p2_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pizero_ptratio_best_vs_pt_barrel_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1Puppi_pizero_ptratio_0p2_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pizero_ptratio_best_vs_pt_endcap_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1Puppi_pizero_ptratio_0p2_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pizero_ptratio_best_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1Puppi_pizero_ptratio_0p2_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pizero_ptratio_best_vs_pt_hf_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          }
          h_L1Puppi_pizero_ptratio_0p2_vs_eta_->Fill(mc_.eta, recopair.second.pt02 / mc_.pt);
          h_L1Puppi_pizero_ptratio_best_vs_eta_->Fill(mc_.eta, recopair.second.ptbest / mc_.pt);
        }
      }
      if (abs(mc_.id) == 211 && fabs(mc_.iso04) < 0.05) {
        if (recopair.first.name() == "L1PF") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1PF_pion_ptratio_0p2_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pion_ptratio_best_vs_pt_barrel_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1PF_pion_ptratio_0p2_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pion_ptratio_best_vs_pt_endcap_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1PF_pion_ptratio_0p2_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pion_ptratio_best_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1PF_pion_ptratio_0p2_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1PF_pion_ptratio_best_vs_pt_hf_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          }
          h_L1PF_pion_ptratio_0p2_vs_eta_->Fill(mc_.eta, recopair.second.pt02 / mc_.pt);
          h_L1PF_pion_ptratio_best_vs_eta_->Fill(mc_.eta, recopair.second.ptbest / mc_.pt);
        }
        if (recopair.first.name() == "L1Puppi") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1Puppi_pion_ptratio_0p2_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pion_ptratio_best_vs_pt_barrel_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1Puppi_pion_ptratio_0p2_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pion_ptratio_best_vs_pt_endcap_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1Puppi_pion_ptratio_0p2_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pion_ptratio_best_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1Puppi_pion_ptratio_0p2_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt02 / mc_.pt);
            h_L1Puppi_pion_ptratio_best_vs_pt_hf_->Fill(mc_.pt, recopair.second.ptbest / mc_.pt);
          }
          h_L1Puppi_pion_ptratio_0p2_vs_eta_->Fill(mc_.eta, recopair.second.pt02 / mc_.pt);
          h_L1Puppi_pion_ptratio_best_vs_eta_->Fill(mc_.eta, recopair.second.ptbest / mc_.pt);
        }
      }
      if (abs(mc_.id) == 0) {
        if (recopair.first.name() == "L1PF") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1PF_jet_ptratio_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1PF_jet_ptratio_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1PF_jet_ptratio_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1PF_jet_ptratio_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          }
          h_L1PF_jet_ptratio_vs_eta_->Fill(mc_.eta, recopair.second.pt / mc_.pt);
        }
        if (recopair.first.name() == "L1Puppi") {
          if (fabs(mc_.eta) < 1.5) {
            h_L1Puppi_jet_ptratio_vs_pt_barrel_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          } else if (fabs(mc_.eta) < 2.5) {
            h_L1Puppi_jet_ptratio_vs_pt_endcap_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          } else if (fabs(mc_.eta) < 3.) {
            h_L1Puppi_jet_ptratio_vs_pt_ecnotk_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          } else if (fabs(mc_.eta) < 5.) {
            h_L1Puppi_jet_ptratio_vs_pt_hf_->Fill(mc_.pt, recopair.second.pt / mc_.pt);
          }
          h_L1Puppi_jet_ptratio_vs_eta_->Fill(mc_.eta, recopair.second.pt / mc_.pt);
        }
      }
    }
  }
  for (auto& recopair : reco_) {
    recopair.first.clear();
    recopair.second.clear();
  }

  for (const auto& pfc : *l1pfs) {
    h_L1PF_pt_->Fill(pfc.pt());
    h_L1PF_eta_->Fill(pfc.eta());
    if (abs(pfc.pdgId()) == 13) {
      h_L1PF_pt_mu_->Fill(pfc.pt());
      h_L1PF_eta_mu_->Fill(pfc.eta());
    } else if (abs(pfc.pdgId()) == 11) {
      h_L1PF_pt_el_->Fill(pfc.pt());
      h_L1PF_eta_el_->Fill(pfc.eta());
    } else if (abs(pfc.pdgId()) == 22) {
      h_L1PF_pt_pho_->Fill(pfc.pt());
      h_L1PF_eta_pho_->Fill(pfc.eta());
    } else if (abs(pfc.pdgId()) == 211) {
      h_L1PF_pt_ch_->Fill(pfc.pt());
      h_L1PF_eta_ch_->Fill(pfc.eta());
    } else if (abs(pfc.pdgId()) == 130) {
      h_L1PF_pt_nh_->Fill(pfc.pt());
      h_L1PF_eta_nh_->Fill(pfc.eta());
    }
  }  // loop over L1 PF
  for (const auto& pupc : *l1pups) {
    h_L1Puppi_pt_->Fill(pupc.pt());
    h_L1Puppi_eta_->Fill(pupc.eta());
    if (abs(pupc.pdgId()) == 13) {
      h_L1Puppi_pt_mu_->Fill(pupc.pt());
      h_L1Puppi_eta_mu_->Fill(pupc.eta());
    } else if (abs(pupc.pdgId()) == 11) {
      h_L1Puppi_pt_el_->Fill(pupc.pt());
      h_L1Puppi_eta_el_->Fill(pupc.eta());
    } else if (abs(pupc.pdgId()) == 22) {
      h_L1Puppi_pt_pho_->Fill(pupc.pt());
      h_L1Puppi_eta_pho_->Fill(pupc.eta());
    } else if (abs(pupc.pdgId()) == 211) {
      h_L1Puppi_pt_ch_->Fill(pupc.pt());
      h_L1Puppi_eta_ch_->Fill(pupc.eta());
    } else if (abs(pupc.pdgId()) == 130) {
      h_L1Puppi_pt_nh_->Fill(pupc.pt());
      h_L1Puppi_eta_nh_->Fill(pupc.eta());
    }
  }  // loop over L1 Puppi
}

//
// -------------------------------------- endRun --------------------------------------------
//
//
// -------------------------------------- book histograms --------------------------------------------
//
void L1TPhase2CorrelatorOffline::bookPhase2CorrelatorHistos(DQMStore::IBooker& ibooker) {
  ibooker.cd();
  ibooker.setCurrentFolder(histFolder_);
  ibooker.setScope(MonitorElementData::Scope::RUN);

  dqmoffline::l1t::HistDefinition resVsPtDef = histDefinitions_[PlotConfig::resVsPt];
  dqmoffline::l1t::HistDefinition resVsEtaDef = histDefinitions_[PlotConfig::resVsEta];
  dqmoffline::l1t::HistDefinition ptDistDef = histDefinitions_[PlotConfig::ptDist];
  dqmoffline::l1t::HistDefinition etaDistDef = histDefinitions_[PlotConfig::etaDist];

  int ptratio_nbins = 300;
  float ptratio_lo = 0.;
  float ptratio_hi = 3.;

  h_L1PF_pt_ = ibooker.book1D("PF_pt", "L1 PF p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1PF_eta_ = ibooker.book1D("PF_eta", "L1 PF #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);
  h_L1Puppi_pt_ = ibooker.book1D("Puppi_pt", "L1 PUPPI p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1Puppi_eta_ = ibooker.book1D("Puppi_eta", "L1 PUPPI #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);

  h_L1PF_pt_mu_ = ibooker.book1D("PF_pt_mu", "L1 PF Muon p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1PF_eta_mu_ = ibooker.book1D("PF_eta_mu", "L1 PF Muon #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);
  h_L1Puppi_pt_mu_ =
      ibooker.book1D("Puppi_pt_mu", "L1 PUPPI Muon p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1Puppi_eta_mu_ =
      ibooker.book1D("Puppi_eta_mu", "L1 PUPPI Muon #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);

  h_L1PF_pt_el_ = ibooker.book1D("PF_pt_el", "L1 PF Electron p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1PF_eta_el_ =
      ibooker.book1D("PF_eta_el", "L1 PF Electron #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);
  h_L1Puppi_pt_el_ =
      ibooker.book1D("Puppi_pt_el", "L1 PUPPI Electron p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1Puppi_eta_el_ =
      ibooker.book1D("Puppi_eta_el", "L1 PUPPI Electron #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);

  h_L1PF_pt_pho_ = ibooker.book1D("PF_pt_pho", "L1 PF Photon p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1PF_eta_pho_ =
      ibooker.book1D("PF_eta_pho", "L1 PF Photon #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);
  h_L1Puppi_pt_pho_ =
      ibooker.book1D("Puppi_pt_pho", "L1 PUPPI Photon p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1Puppi_eta_pho_ =
      ibooker.book1D("Puppi_eta_pho", "L1 PUPPI Photon #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);

  h_L1PF_pt_ch_ =
      ibooker.book1D("PF_pt_ch", "L1 PF Charged Hadron p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1PF_eta_ch_ =
      ibooker.book1D("PF_eta_ch", "L1 PF Charged Hadron #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);
  h_L1Puppi_pt_ch_ =
      ibooker.book1D("Puppi_pt_ch", "L1 PUPPI Charged Hadron p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1Puppi_eta_ch_ = ibooker.book1D(
      "Puppi_eta_ch", "L1 PUPPI Charged Hadron #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);

  h_L1PF_pt_nh_ =
      ibooker.book1D("PF_pt_nh", "L1 PF Neutral Hadron p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1PF_eta_nh_ =
      ibooker.book1D("PF_eta_nh", "L1 PF Neutral Hadron #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);
  h_L1Puppi_pt_nh_ =
      ibooker.book1D("Puppi_pt_nh", "L1 PUPPI Neutral Hadron p_{T}", ptDistDef.nbinsX, ptDistDef.xmin, ptDistDef.xmax);
  h_L1Puppi_eta_nh_ = ibooker.book1D(
      "Puppi_eta_nh", "L1 PUPPI Neutral Hadron #eta", etaDistDef.nbinsX, etaDistDef.xmin, etaDistDef.xmax);

  ibooker.setCurrentFolder(respresolFolder_);

  h_L1PF_electron_ptratio_0p2_vs_pt_barrel_ = ibooker.book2D("L1PFElectronPtRatio0p2VsPtBarrel",
                                                             "L1 PF Electron L1/Gen (#Delta R < 0.2) vs p_{T}, Barrel",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax,
                                                             ptratio_nbins,
                                                             ptratio_lo,
                                                             ptratio_hi);

  h_L1PF_electron_ptratio_0p2_vs_pt_endcap_ = ibooker.book2D("L1PFElectronPtRatio0p2VsPtEndcap",
                                                             "L1 PF Electron L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax,
                                                             ptratio_nbins,
                                                             ptratio_lo,
                                                             ptratio_hi);

  h_L1PF_electron_ptratio_0p2_vs_pt_ecnotk_ =
      ibooker.book2D("L1PFElectronPtRatio0p2VsPtEndcapNoTk",
                     "L1 PF Electron L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1PF_electron_ptratio_0p2_vs_pt_hf_ = ibooker.book2D("L1PFElectronPtRatio0p2VsPtHF",
                                                         "L1 PF Electron L1/Gen (#Delta R < 0.2) vs p_{T}, HF",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax,
                                                         ptratio_nbins,
                                                         ptratio_lo,
                                                         ptratio_hi);

  h_L1PF_electron_ptratio_0p2_vs_eta_ = ibooker.book2D("L1PFElectronPtRatio0p2VsEta",
                                                       "L1 PF Electron L1/Gen (#Delta R < 0.2) vs #eta",
                                                       resVsEtaDef.nbinsX,
                                                       resVsEtaDef.xmin,
                                                       resVsEtaDef.xmax,
                                                       ptratio_nbins,
                                                       ptratio_lo,
                                                       ptratio_hi);

  h_L1Puppi_electron_ptratio_0p2_vs_pt_barrel_ =
      ibooker.book2D("L1PUPPIElectronPtRatio0p2VsPtBarrel",
                     "L1 PUPPI Electron L1/Gen (#Delta R < 0.2) vs p_{T}, Barrel",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1Puppi_electron_ptratio_0p2_vs_pt_endcap_ =
      ibooker.book2D("L1PUPPIElectronPtRatio0p2VsPtEndcap",
                     "L1 PUPPI Electron L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1Puppi_electron_ptratio_0p2_vs_pt_ecnotk_ =
      ibooker.book2D("L1PUPPIElectronPtRatio0p2VsPtEndcapNoTk",
                     "L1 PUPPI Electron L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1Puppi_electron_ptratio_0p2_vs_pt_hf_ = ibooker.book2D("L1PUPPIElectronPtRatio0p2VsPtHF",
                                                            "L1 PUPPI Electron L1/Gen (#Delta R < 0.2) vs p_{T}, HF",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax,
                                                            ptratio_nbins,
                                                            ptratio_lo,
                                                            ptratio_hi);

  h_L1Puppi_electron_ptratio_0p2_vs_eta_ = ibooker.book2D("L1PUPPIElectronPtRatio0p2VsEta",
                                                          "L1 PUPPI Electron L1/Gen (#Delta R < 0.2) vs #eta",
                                                          resVsEtaDef.nbinsX,
                                                          resVsEtaDef.xmin,
                                                          resVsEtaDef.xmax,
                                                          ptratio_nbins,
                                                          ptratio_lo,
                                                          ptratio_hi);

  h_L1PF_electron_ptratio_best_vs_pt_barrel_ = ibooker.book2D("L1PFElectronPtRatioBestVsPtBarrel",
                                                              "L1 PF Electron L1/Gen (Best) vs p_{T}, Barrel",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax,
                                                              ptratio_nbins,
                                                              ptratio_lo,
                                                              ptratio_hi);

  h_L1PF_electron_ptratio_best_vs_pt_endcap_ = ibooker.book2D("L1PFElectronPtRatioBestVsPtEndcap",
                                                              "L1 PF Electron L1/Gen (Best) vs p_{T}, Endcap",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax,
                                                              ptratio_nbins,
                                                              ptratio_lo,
                                                              ptratio_hi);

  h_L1PF_electron_ptratio_best_vs_pt_ecnotk_ = ibooker.book2D("L1PFElectronPtRatioBestVsPtEndcapNoTk",
                                                              "L1 PF Electron L1/Gen (Best) vs p_{T}, Endcap No Tk",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax,
                                                              ptratio_nbins,
                                                              ptratio_lo,
                                                              ptratio_hi);

  h_L1PF_electron_ptratio_best_vs_pt_hf_ = ibooker.book2D("L1PFElectronPtRatioBestVsPtHF",
                                                          "L1 PF Electron L1/Gen (Best) vs p_{T}, HF",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax,
                                                          ptratio_nbins,
                                                          ptratio_lo,
                                                          ptratio_hi);

  h_L1PF_electron_ptratio_best_vs_eta_ = ibooker.book2D("L1PFElectronPtRatioBestVsEta",
                                                        "L1 PF Electron L1/Gen (Best) vs #eta",
                                                        resVsEtaDef.nbinsX,
                                                        resVsEtaDef.xmin,
                                                        resVsEtaDef.xmax,
                                                        ptratio_nbins,
                                                        ptratio_lo,
                                                        ptratio_hi);

  h_L1Puppi_electron_ptratio_best_vs_pt_barrel_ = ibooker.book2D("L1PUPPIElectronPtRatioBestVsPtBarrel",
                                                                 "L1 PUPPI Electron L1/Gen (Best) vs p_{T}, Barrel",
                                                                 resVsPtDef.nbinsX,
                                                                 resVsPtDef.xmin,
                                                                 resVsPtDef.xmax,
                                                                 ptratio_nbins,
                                                                 ptratio_lo,
                                                                 ptratio_hi);

  h_L1Puppi_electron_ptratio_best_vs_pt_endcap_ = ibooker.book2D("L1PUPPIElectronPtRatioBestVsPtEndcap",
                                                                 "L1 PUPPI Electron L1/Gen (Best) vs p_{T}, Endcap",
                                                                 resVsPtDef.nbinsX,
                                                                 resVsPtDef.xmin,
                                                                 resVsPtDef.xmax,
                                                                 ptratio_nbins,
                                                                 ptratio_lo,
                                                                 ptratio_hi);

  h_L1Puppi_electron_ptratio_best_vs_pt_ecnotk_ =
      ibooker.book2D("L1PUPPIElectronPtRatioBestVsPtEndcapNoTk",
                     "L1 PUPPI Electron L1/Gen (Best) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1Puppi_electron_ptratio_best_vs_pt_hf_ = ibooker.book2D("L1PUPPIElectronPtRatioBestVsPtHF",
                                                             "L1 PUPPI Electron L1/Gen (Best) vs p_{T}, HF",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax,
                                                             ptratio_nbins,
                                                             ptratio_lo,
                                                             ptratio_hi);

  h_L1Puppi_electron_ptratio_best_vs_eta_ = ibooker.book2D("L1PUPPIElectronPtRatioBestVsEta",
                                                           "L1 PUPPI Electron L1/Gen (Best) vs #eta",
                                                           resVsEtaDef.nbinsX,
                                                           resVsEtaDef.xmin,
                                                           resVsEtaDef.xmax,
                                                           ptratio_nbins,
                                                           ptratio_lo,
                                                           ptratio_hi);

  h_L1PF_pizero_ptratio_0p2_vs_pt_barrel_ = ibooker.book2D("L1PFPi0PtRatio0p2VsPtBarrel",
                                                           "L1 PF Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, Barrel",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax,
                                                           ptratio_nbins,
                                                           ptratio_lo,
                                                           ptratio_hi);

  h_L1PF_pizero_ptratio_0p2_vs_pt_endcap_ = ibooker.book2D("L1PFPi0PtRatio0p2VsPtEndcap",
                                                           "L1 PF Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax,
                                                           ptratio_nbins,
                                                           ptratio_lo,
                                                           ptratio_hi);

  h_L1PF_pizero_ptratio_0p2_vs_pt_ecnotk_ = ibooker.book2D("L1PFPi0PtRatio0p2VsPtEndcapNoTk",
                                                           "L1 PF Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax,
                                                           ptratio_nbins,
                                                           ptratio_lo,
                                                           ptratio_hi);

  h_L1PF_pizero_ptratio_0p2_vs_pt_hf_ = ibooker.book2D("L1PFPi0PtRatio0p2VsPtHF",
                                                       "L1 PF Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, HF",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax,
                                                       ptratio_nbins,
                                                       ptratio_lo,
                                                       ptratio_hi);

  h_L1PF_pizero_ptratio_0p2_vs_eta_ = ibooker.book2D("L1PFPi0PtRatio0p2VsEta",
                                                     "L1 PF Pi0 L1/Gen (#Delta R < 0.2) vs #eta",
                                                     resVsEtaDef.nbinsX,
                                                     resVsEtaDef.xmin,
                                                     resVsEtaDef.xmax,
                                                     ptratio_nbins,
                                                     ptratio_lo,
                                                     ptratio_hi);

  h_L1Puppi_pizero_ptratio_0p2_vs_pt_barrel_ = ibooker.book2D("L1PUPPIPi0PtRatio0p2VsPtBarrel",
                                                              "L1 PUPPI Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, Barrel",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax,
                                                              ptratio_nbins,
                                                              ptratio_lo,
                                                              ptratio_hi);

  h_L1Puppi_pizero_ptratio_0p2_vs_pt_endcap_ = ibooker.book2D("L1PUPPIPi0PtRatio0p2VsPtEndcap",
                                                              "L1 PUPPI Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax,
                                                              ptratio_nbins,
                                                              ptratio_lo,
                                                              ptratio_hi);

  h_L1Puppi_pizero_ptratio_0p2_vs_pt_ecnotk_ =
      ibooker.book2D("L1PUPPIPi0PtRatio0p2VsPtEndcapNoTk",
                     "L1 PUPPI Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1Puppi_pizero_ptratio_0p2_vs_pt_hf_ = ibooker.book2D("L1PUPPIPi0PtRatio0p2VsPtHF",
                                                          "L1 PUPPI Pi0 L1/Gen (#Delta R < 0.2) vs p_{T}, HF",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax,
                                                          ptratio_nbins,
                                                          ptratio_lo,
                                                          ptratio_hi);

  h_L1Puppi_pizero_ptratio_0p2_vs_eta_ = ibooker.book2D("L1PUPPIPi0PtRatio0p2VsEta",
                                                        "L1 PUPPI Pi0 L1/Gen (#Delta R < 0.2) vs #eta",
                                                        resVsEtaDef.nbinsX,
                                                        resVsEtaDef.xmin,
                                                        resVsEtaDef.xmax,
                                                        ptratio_nbins,
                                                        ptratio_lo,
                                                        ptratio_hi);

  h_L1PF_pizero_ptratio_best_vs_pt_barrel_ = ibooker.book2D("L1PFPi0PtRatioBestVsPtBarrel",
                                                            "L1 PF Pi0 L1/Gen (Best) vs p_{T}, Barrel",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax,
                                                            ptratio_nbins,
                                                            ptratio_lo,
                                                            ptratio_hi);

  h_L1PF_pizero_ptratio_best_vs_pt_endcap_ = ibooker.book2D("L1PFPi0PtRatioBestVsPtEndcap",
                                                            "L1 PF Pi0 L1/Gen (Best) vs p_{T}, Endcap",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax,
                                                            ptratio_nbins,
                                                            ptratio_lo,
                                                            ptratio_hi);

  h_L1PF_pizero_ptratio_best_vs_pt_ecnotk_ = ibooker.book2D("L1PFPi0PtRatioBestVsPtEndcapNoTk",
                                                            "L1 PF Pi0 L1/Gen (Best) vs p_{T}, Endcap No Tk",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax,
                                                            ptratio_nbins,
                                                            ptratio_lo,
                                                            ptratio_hi);

  h_L1PF_pizero_ptratio_best_vs_pt_hf_ = ibooker.book2D("L1PFPi0PtRatioBestVsPtHF",
                                                        "L1 PF Pi0 L1/Gen (Best) vs p_{T}, HF",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax,
                                                        ptratio_nbins,
                                                        ptratio_lo,
                                                        ptratio_hi);

  h_L1PF_pizero_ptratio_best_vs_eta_ = ibooker.book2D("L1PFPi0PtRatioBestVsEta",
                                                      "L1 PF Pi0 L1/Gen (Best) vs #eta",
                                                      resVsEtaDef.nbinsX,
                                                      resVsEtaDef.xmin,
                                                      resVsEtaDef.xmax,
                                                      ptratio_nbins,
                                                      ptratio_lo,
                                                      ptratio_hi);

  h_L1Puppi_pizero_ptratio_best_vs_pt_barrel_ = ibooker.book2D("L1PUPPIPi0PtRatioBestVsPtBarrel",
                                                               "L1 PUPPI Pi0 L1/Gen (Best) vs p_{T}, Barrel",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax,
                                                               ptratio_nbins,
                                                               ptratio_lo,
                                                               ptratio_hi);

  h_L1Puppi_pizero_ptratio_best_vs_pt_endcap_ = ibooker.book2D("L1PUPPIPi0PtRatioBestVsPtEndcap",
                                                               "L1 PUPPI Pi0 L1/Gen (Best) vs p_{T}, Endcap",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax,
                                                               ptratio_nbins,
                                                               ptratio_lo,
                                                               ptratio_hi);

  h_L1Puppi_pizero_ptratio_best_vs_pt_ecnotk_ = ibooker.book2D("L1PUPPIPi0PtRatioBestVsPtEndcapNoTk",
                                                               "L1 PUPPI Pi0 L1/Gen (Best) vs p_{T}, Endcap No Tk",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax,
                                                               ptratio_nbins,
                                                               ptratio_lo,
                                                               ptratio_hi);

  h_L1Puppi_pizero_ptratio_best_vs_pt_hf_ = ibooker.book2D("L1PUPPIPi0PtRatioBestVsPtHF",
                                                           "L1 PUPPI Pi0 L1/Gen (Best) vs p_{T}, HF",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax,
                                                           ptratio_nbins,
                                                           ptratio_lo,
                                                           ptratio_hi);

  h_L1Puppi_pizero_ptratio_best_vs_eta_ = ibooker.book2D("L1PUPPIPi0PtRatioBestVsEta",
                                                         "L1 PUPPI Pi0 L1/Gen (Best) vs #eta",
                                                         resVsEtaDef.nbinsX,
                                                         resVsEtaDef.xmin,
                                                         resVsEtaDef.xmax,
                                                         ptratio_nbins,
                                                         ptratio_lo,
                                                         ptratio_hi);

  h_L1PF_pion_ptratio_0p2_vs_pt_barrel_ = ibooker.book2D("L1PFPionPtRatio0p2VsPtBarrel",
                                                         "L1 PF Pion L1/Gen (#Delta R < 0.2) vs p_{T}, Barrel",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax,
                                                         ptratio_nbins,
                                                         ptratio_lo,
                                                         ptratio_hi);

  h_L1PF_pion_ptratio_0p2_vs_pt_endcap_ = ibooker.book2D("L1PFPionPtRatio0p2VsPtEndcap",
                                                         "L1 PF Pion L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax,
                                                         ptratio_nbins,
                                                         ptratio_lo,
                                                         ptratio_hi);

  h_L1PF_pion_ptratio_0p2_vs_pt_ecnotk_ = ibooker.book2D("L1PFPionPtRatio0p2VsPtEndcapNoTk",
                                                         "L1 PF Pion L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax,
                                                         ptratio_nbins,
                                                         ptratio_lo,
                                                         ptratio_hi);

  h_L1PF_pion_ptratio_0p2_vs_pt_hf_ = ibooker.book2D("L1PFPionPtRatio0p2VsPtHF",
                                                     "L1 PF Pion L1/Gen (#Delta R < 0.2) vs p_{T}, HF",
                                                     resVsPtDef.nbinsX,
                                                     resVsPtDef.xmin,
                                                     resVsPtDef.xmax,
                                                     ptratio_nbins,
                                                     ptratio_lo,
                                                     ptratio_hi);

  h_L1PF_pion_ptratio_0p2_vs_eta_ = ibooker.book2D("L1PFPionPtRatio0p2VsEta",
                                                   "L1 PF Pion L1/Gen (#Delta R < 0.2) vs #eta",
                                                   resVsEtaDef.nbinsX,
                                                   resVsEtaDef.xmin,
                                                   resVsEtaDef.xmax,
                                                   ptratio_nbins,
                                                   ptratio_lo,
                                                   ptratio_hi);

  h_L1Puppi_pion_ptratio_0p2_vs_pt_barrel_ = ibooker.book2D("L1PUPPIPionPtRatio0p2VsPtBarrel",
                                                            "L1 PUPPI Pion L1/Gen (#Delta R < 0.2) vs p_{T}, Barrel",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax,
                                                            ptratio_nbins,
                                                            ptratio_lo,
                                                            ptratio_hi);

  h_L1Puppi_pion_ptratio_0p2_vs_pt_endcap_ = ibooker.book2D("L1PUPPIPionPtRatio0p2VsPtEndcap",
                                                            "L1 PUPPI Pion L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax,
                                                            ptratio_nbins,
                                                            ptratio_lo,
                                                            ptratio_hi);

  h_L1Puppi_pion_ptratio_0p2_vs_pt_ecnotk_ =
      ibooker.book2D("L1PUPPIPionPtRatio0p2VsPtEndcapNoTk",
                     "L1 PUPPI Pion L1/Gen (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax,
                     ptratio_nbins,
                     ptratio_lo,
                     ptratio_hi);

  h_L1Puppi_pion_ptratio_0p2_vs_pt_hf_ = ibooker.book2D("L1PUPPIPionPtRatio0p2VsPtHF",
                                                        "L1 PUPPI Pion L1/Gen (#Delta R < 0.2) vs p_{T}, HF",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax,
                                                        ptratio_nbins,
                                                        ptratio_lo,
                                                        ptratio_hi);

  h_L1Puppi_pion_ptratio_0p2_vs_eta_ = ibooker.book2D("L1PUPPIPionPtRatio0p2VsEta",
                                                      "L1 PUPPI Pion L1/Gen (#Delta R < 0.2) vs #eta",
                                                      resVsEtaDef.nbinsX,
                                                      resVsEtaDef.xmin,
                                                      resVsEtaDef.xmax,
                                                      ptratio_nbins,
                                                      ptratio_lo,
                                                      ptratio_hi);

  h_L1PF_pion_ptratio_best_vs_pt_barrel_ = ibooker.book2D("L1PFPionPtRatioBestVsPtBarrel",
                                                          "L1 PF Pion L1/Gen (Best) vs p_{T}, Barrel",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax,
                                                          ptratio_nbins,
                                                          ptratio_lo,
                                                          ptratio_hi);

  h_L1PF_pion_ptratio_best_vs_pt_endcap_ = ibooker.book2D("L1PFPionPtRatioBestVsPtEndcap",
                                                          "L1 PF Pion L1/Gen (Best) vs p_{T}, Endcap",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax,
                                                          ptratio_nbins,
                                                          ptratio_lo,
                                                          ptratio_hi);

  h_L1PF_pion_ptratio_best_vs_pt_ecnotk_ = ibooker.book2D("L1PFPionPtRatioBestVsPtEndcapNoTk",
                                                          "L1 PF Pion L1/Gen (Best) vs p_{T}, Endcap No Tk",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax,
                                                          ptratio_nbins,
                                                          ptratio_lo,
                                                          ptratio_hi);

  h_L1PF_pion_ptratio_best_vs_pt_hf_ = ibooker.book2D("L1PFPionPtRatioBestVsPtHF",
                                                      "L1 PF Pion L1/Gen (Best) vs p_{T}, HF",
                                                      resVsPtDef.nbinsX,
                                                      resVsPtDef.xmin,
                                                      resVsPtDef.xmax,
                                                      ptratio_nbins,
                                                      ptratio_lo,
                                                      ptratio_hi);

  h_L1PF_pion_ptratio_best_vs_eta_ = ibooker.book2D("L1PFPionPtRatioBestVsEta",
                                                    "L1 PF Pion L1/Gen (Best) vs #eta",
                                                    resVsEtaDef.nbinsX,
                                                    resVsEtaDef.xmin,
                                                    resVsEtaDef.xmax,
                                                    ptratio_nbins,
                                                    ptratio_lo,
                                                    ptratio_hi);

  h_L1Puppi_pion_ptratio_best_vs_pt_barrel_ = ibooker.book2D("L1PUPPIPionPtRatioBestVsPtBarrel",
                                                             "L1 PUPPI Pion L1/Gen (Best) vs p_{T}, Barrel",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax,
                                                             ptratio_nbins,
                                                             ptratio_lo,
                                                             ptratio_hi);

  h_L1Puppi_pion_ptratio_best_vs_pt_endcap_ = ibooker.book2D("L1PUPPIPionPtRatioBestVsPtEndcap",
                                                             "L1 PUPPI Pion L1/Gen (Best) vs p_{T}, Endcap",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax,
                                                             ptratio_nbins,
                                                             ptratio_lo,
                                                             ptratio_hi);

  h_L1Puppi_pion_ptratio_best_vs_pt_ecnotk_ = ibooker.book2D("L1PUPPIPionPtRatioBestVsPtEndcapNoTk",
                                                             "L1 PUPPI Pion L1/Gen (Best) vs p_{T}, Endcap No Tk",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax,
                                                             ptratio_nbins,
                                                             ptratio_lo,
                                                             ptratio_hi);

  h_L1Puppi_pion_ptratio_best_vs_pt_hf_ = ibooker.book2D("L1PUPPIPionPtRatioBestVsPtHF",
                                                         "L1 PUPPI Pion L1/Gen (Best) vs p_{T}, HF",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax,
                                                         ptratio_nbins,
                                                         ptratio_lo,
                                                         ptratio_hi);

  h_L1Puppi_pion_ptratio_best_vs_eta_ = ibooker.book2D("L1PUPPIPionPtRatioBestVsEta",
                                                       "L1 PUPPI Pion L1/Gen (Best) vs #eta",
                                                       resVsEtaDef.nbinsX,
                                                       resVsEtaDef.xmin,
                                                       resVsEtaDef.xmax,
                                                       ptratio_nbins,
                                                       ptratio_lo,
                                                       ptratio_hi);

  h_L1PF_jet_ptratio_vs_pt_barrel_ = ibooker.book2D("L1PFJetPtRatioVsPtBarrel",
                                                    "L1 PF Jet L1/Gen vs p_{T}, Barrel",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax,
                                                    ptratio_nbins,
                                                    ptratio_lo,
                                                    ptratio_hi);

  h_L1PF_jet_ptratio_vs_pt_endcap_ = ibooker.book2D("L1PFJetPtRatioVsPtEndcap",
                                                    "L1 PF Jet L1/Gen vs p_{T}, Endcap",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax,
                                                    ptratio_nbins,
                                                    ptratio_lo,
                                                    ptratio_hi);

  h_L1PF_jet_ptratio_vs_pt_ecnotk_ = ibooker.book2D("L1PFJetPtRatioVsPtEndcapNoTk",
                                                    "L1 PF Jet L1/Gen vs p_{T}, Endcap No Tk",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax,
                                                    ptratio_nbins,
                                                    ptratio_lo,
                                                    ptratio_hi);

  h_L1PF_jet_ptratio_vs_pt_hf_ = ibooker.book2D("L1PFJetPtRatioVsPtHF",
                                                "L1 PF Jet L1/Gen vs p_{T}, HF",
                                                resVsPtDef.nbinsX,
                                                resVsPtDef.xmin,
                                                resVsPtDef.xmax,
                                                ptratio_nbins,
                                                ptratio_lo,
                                                ptratio_hi);

  h_L1PF_jet_ptratio_vs_eta_ = ibooker.book2D("L1PFJetPtRatioVsEta",
                                              "L1 PF Jet L1/Gen vs #eta",
                                              resVsEtaDef.nbinsX,
                                              resVsEtaDef.xmin,
                                              resVsEtaDef.xmax,
                                              ptratio_nbins,
                                              ptratio_lo,
                                              ptratio_hi);

  h_L1Puppi_jet_ptratio_vs_pt_barrel_ = ibooker.book2D("L1PUPPIJetPtRatioVsPtBarrel",
                                                       "L1 PUPPI Jet L1/Gen vs p_{T}, Barrel",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax,
                                                       ptratio_nbins,
                                                       ptratio_lo,
                                                       ptratio_hi);

  h_L1Puppi_jet_ptratio_vs_pt_endcap_ = ibooker.book2D("L1PUPPIJetPtRatioVsPtEndcap",
                                                       "L1 PUPPI Jet L1/Gen vs p_{T}, Endcap",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax,
                                                       ptratio_nbins,
                                                       ptratio_lo,
                                                       ptratio_hi);

  h_L1Puppi_jet_ptratio_vs_pt_ecnotk_ = ibooker.book2D("L1PUPPIJetPtRatioVsPtEndcapNoTk",
                                                       "L1 PUPPI Jet L1/Gen vs p_{T}, EndcapNoTk",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax,
                                                       ptratio_nbins,
                                                       ptratio_lo,
                                                       ptratio_hi);

  h_L1Puppi_jet_ptratio_vs_pt_hf_ = ibooker.book2D("L1PUPPIJetPtRatioVsPtHF",
                                                   "L1 PUPPI Jet L1/Gen vs p_{T}, HF",
                                                   resVsPtDef.nbinsX,
                                                   resVsPtDef.xmin,
                                                   resVsPtDef.xmax,
                                                   ptratio_nbins,
                                                   ptratio_lo,
                                                   ptratio_hi);

  h_L1Puppi_jet_ptratio_vs_eta_ = ibooker.book2D("L1PUPPIJetPtRatioVsEta",
                                                 "L1 PUPPI Jet L1/Gen vs #eta",
                                                 resVsEtaDef.nbinsX,
                                                 resVsEtaDef.xmin,
                                                 resVsEtaDef.xmax,
                                                 ptratio_nbins,
                                                 ptratio_lo,
                                                 ptratio_hi);

  ibooker.setCurrentFolder(histFolder_);
  //Response
  h_L1PF_electron_response_0p2_pt_barrel_ = ibooker.book1D("L1PFElectronResponse0p2VsPtBarrel",
                                                           "L1 PF Electron Response (#Delta R < 0.2) vs p_{T}, Barrel",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1PF_electron_response_0p2_pt_endcap_ = ibooker.book1D("L1PFElectronResponse0p2VsPtEndcap",
                                                           "L1 PF Electron Response (#Delta R < 0.2) vs p_{T}, Endcap",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1PF_electron_response_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PFElectronResponse0p2VsPtEndcapNoTk",
                     "L1 PF Electron Response (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_electron_response_0p2_pt_hf_ = ibooker.book1D("L1PFElectronResponse0p2VsPtHF",
                                                       "L1 PF Electron Response (#Delta R < 0.2) vs p_{T}, HF",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1PF_electron_response_0p2_eta_ = ibooker.book1D("L1PFElectronResponse0p2VsEta",
                                                     "L1 PF Electron Response (#Delta R < 0.2) vs #eta",
                                                     resVsEtaDef.nbinsX,
                                                     resVsEtaDef.xmin,
                                                     resVsEtaDef.xmax);

  h_L1Puppi_electron_response_0p2_pt_barrel_ =
      ibooker.book1D("L1PUPPIElectronResponse0p2VsPtBarrel",
                     "L1 PUPPI Electron Response (#Delta R < 0.2) vs p_{T}, Barrel",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_response_0p2_pt_endcap_ =
      ibooker.book1D("L1PUPPIElectronResponse0p2VsPtEndcap",
                     "L1 PUPPI Electron Response (#Delta R < 0.2) vs p_{T}, Endcap",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_response_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIElectronResponse0p2VsPtEndcapNoTk",
                     "L1 PUPPI Electron Response (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_response_0p2_pt_hf_ = ibooker.book1D("L1PUPPIElectronResponse0p2VsPtHF",
                                                          "L1 PUPPI Electron Response (#Delta R < 0.2) vs p_{T}, HF",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1Puppi_electron_response_0p2_eta_ = ibooker.book1D("L1PUPPIElectronResponse0p2VsEta",
                                                        "L1 PUPPI Electron Response (#Delta R < 0.2) vs #eta",
                                                        resVsEtaDef.nbinsX,
                                                        resVsEtaDef.xmin,
                                                        resVsEtaDef.xmax);

  h_L1PF_electron_response_best_pt_barrel_ = ibooker.book1D("L1PFElectronResponseBestVsPtBarrel",
                                                            "L1 PF Electron Response (Best) vs p_{T}, Barrel",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1PF_electron_response_best_pt_endcap_ = ibooker.book1D("L1PFElectronResponseBestVsPtEndcap",
                                                            "L1 PF Electron Response (Best) vs p_{T}, Endcap",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1PF_electron_response_best_pt_ecnotk_ = ibooker.book1D("L1PFElectronResponseBestVsPtEndcapNoTk",
                                                            "L1 PF Electron Response (Best) vs p_{T}, Endcap No Tk",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1PF_electron_response_best_pt_hf_ = ibooker.book1D("L1PFElectronResponseBestVsPtHF",
                                                        "L1 PF Electron Response (Best) vs p_{T}, HF",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1PF_electron_response_best_eta_ = ibooker.book1D("L1PFElectronResponseBestVsEta",
                                                      "L1 PF Electron Response (Best) vs #eta",
                                                      resVsEtaDef.nbinsX,
                                                      resVsEtaDef.xmin,
                                                      resVsEtaDef.xmax);

  h_L1Puppi_electron_response_best_pt_barrel_ = ibooker.book1D("L1PUPPIElectronResponseBestVsPtBarrel",
                                                               "L1 PUPPI Electron Response (Best) vs p_{T}, Barrel",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax);

  h_L1Puppi_electron_response_best_pt_endcap_ = ibooker.book1D("L1PUPPIElectronResponseBestVsPtEndcap",
                                                               "L1 PUPPI Electron Response (Best) vs p_{T}, Endcap",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax);

  h_L1Puppi_electron_response_best_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIElectronResponseBestVsPtEndcapNoTk",
                     "L1 PUPPI Electron Response (Best) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_response_best_pt_hf_ = ibooker.book1D("L1PUPPIElectronResponseBestVsPtHF",
                                                           "L1 PUPPI Electron Response (Best) vs p_{T}, HF",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1Puppi_electron_response_best_eta_ = ibooker.book1D("L1PUPPIElectronResponseBestVsEta",
                                                         "L1 PUPPI Electron Response (Best) vs #eta",
                                                         resVsEtaDef.nbinsX,
                                                         resVsEtaDef.xmin,
                                                         resVsEtaDef.xmax);

  h_L1PF_pizero_response_0p2_pt_barrel_ = ibooker.book1D("L1PFPi0Response0p2VsPtBarrel",
                                                         "L1 PF Pi0 Response (#Delta R < 0.2) vs p_{T}, Barrel",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1PF_pizero_response_0p2_pt_endcap_ = ibooker.book1D("L1PFPi0Response0p2VsPtEndcap",
                                                         "L1 PF Pi0 Response (#Delta R < 0.2) vs p_{T}, Endcap",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1PF_pizero_response_0p2_pt_ecnotk_ = ibooker.book1D("L1PFPi0Response0p2VsPtEndcapNoTk",
                                                         "L1 PF Pi0 Response (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1PF_pizero_response_0p2_pt_hf_ = ibooker.book1D("L1PFPi0Response0p2VsPtHF",
                                                     "L1 PF Pi0 Response (#Delta R < 0.2) vs p_{T}, HF",
                                                     resVsPtDef.nbinsX,
                                                     resVsPtDef.xmin,
                                                     resVsPtDef.xmax);

  h_L1PF_pizero_response_0p2_eta_ = ibooker.book1D("L1PFPi0Response0p2VsEta",
                                                   "L1 PF Pi0 Response (#Delta R < 0.2) vs #eta",
                                                   resVsEtaDef.nbinsX,
                                                   resVsEtaDef.xmin,
                                                   resVsEtaDef.xmax);

  h_L1Puppi_pizero_response_0p2_pt_barrel_ = ibooker.book1D("L1PUPPIPi0Response0p2VsPtBarrel",
                                                            "L1 PUPPI Pi0 Response (#Delta R < 0.2) vs p_{T}, Barrel",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1Puppi_pizero_response_0p2_pt_endcap_ = ibooker.book1D("L1PUPPIPi0Response0p2VsPtEndcap",
                                                            "L1 PUPPI Pi0 Response (#Delta R < 0.2) vs p_{T}, Endcap",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1Puppi_pizero_response_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIPi0Response0p2VsPtEndcapNoTk",
                     "L1 PUPPI Pi0 Response (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pizero_response_0p2_pt_hf_ = ibooker.book1D("L1PUPPIPi0Response0p2VsPtHF",
                                                        "L1 PUPPI Pi0 Response (#Delta R < 0.2) vs p_{T}, HF",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1Puppi_pizero_response_0p2_eta_ = ibooker.book1D("L1PUPPIPi0Response0p2VsEta",
                                                      "L1 PUPPI Pi0 Response (#Delta R < 0.2) vs #eta",
                                                      resVsEtaDef.nbinsX,
                                                      resVsEtaDef.xmin,
                                                      resVsEtaDef.xmax);

  h_L1PF_pizero_response_best_pt_barrel_ = ibooker.book1D("L1PFPi0ResponseBestVsPtBarrel",
                                                          "L1 PF Pi0 Response (Best) vs p_{T}, Barrel",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pizero_response_best_pt_endcap_ = ibooker.book1D("L1PFPi0ResponseBestVsPtEndcap",
                                                          "L1 PF Pi0 Response (Best) vs p_{T}, Endcap",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pizero_response_best_pt_ecnotk_ = ibooker.book1D("L1PFPi0ResponseBestVsPtEndcapNoTk",
                                                          "L1 PF Pi0 Response (Best) vs p_{T}, Endcap No Tk",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pizero_response_best_pt_hf_ = ibooker.book1D("L1PFPi0ResponseBestVsPtHF",
                                                      "L1 PF Pi0 Response (Best) vs p_{T}, HF",
                                                      resVsPtDef.nbinsX,
                                                      resVsPtDef.xmin,
                                                      resVsPtDef.xmax);

  h_L1PF_pizero_response_best_eta_ = ibooker.book1D("L1PFPi0ResponseBestVsEta",
                                                    "L1 PF Pi0 Response (Best) vs #eta",
                                                    resVsEtaDef.nbinsX,
                                                    resVsEtaDef.xmin,
                                                    resVsEtaDef.xmax);

  h_L1Puppi_pizero_response_best_pt_barrel_ = ibooker.book1D("L1PUPPIPi0ResponseBestVsPtBarrel",
                                                             "L1 PUPPI Pi0 Response (Best) vs p_{T}, Barrel",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1Puppi_pizero_response_best_pt_endcap_ = ibooker.book1D("L1PUPPIPi0ResponseBestVsPtEndcap",
                                                             "L1 PUPPI Pi0 Response (Best) vs p_{T}, Endcap",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1Puppi_pizero_response_best_pt_ecnotk_ = ibooker.book1D("L1PUPPIPi0ResponseBestVsPtEndcapNoTk",
                                                             "L1 PUPPI Pi0 Response (Best) vs p_{T}, Endcap No Tk",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1Puppi_pizero_response_best_pt_hf_ = ibooker.book1D("L1PUPPIPi0ResponseBestVsPtHF",
                                                         "L1 PUPPI Pi0 Response (Best) vs p_{T}, HF",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1Puppi_pizero_response_best_eta_ = ibooker.book1D("L1PUPPIPi0ResponseBestVsEta",
                                                       "L1 PUPPI Pi0 Response (Best) vs #eta",
                                                       resVsEtaDef.nbinsX,
                                                       resVsEtaDef.xmin,
                                                       resVsEtaDef.xmax);

  h_L1PF_pion_response_0p2_pt_barrel_ = ibooker.book1D("L1PFPionResponse0p2VsPtBarrel",
                                                       "L1 PF Pion Response (#Delta R < 0.2) vs p_{T}, Barrel",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1PF_pion_response_0p2_pt_endcap_ = ibooker.book1D("L1PFPionResponse0p2VsPtEndcap",
                                                       "L1 PF Pion Response (#Delta R < 0.2) vs p_{T}, Endcap",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1PF_pion_response_0p2_pt_ecnotk_ = ibooker.book1D("L1PFPionResponse0p2VsPtEndcapNoTk",
                                                       "L1 PF Pion Response (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1PF_pion_response_0p2_pt_hf_ = ibooker.book1D("L1PFPionResponse0p2VsPtHF",
                                                   "L1 PF Pion Response (#Delta R < 0.2) vs p_{T}, HF",
                                                   resVsPtDef.nbinsX,
                                                   resVsPtDef.xmin,
                                                   resVsPtDef.xmax);

  h_L1PF_pion_response_0p2_eta_ = ibooker.book1D("L1PFPionResponse0p2VsEta",
                                                 "L1 PF Pion Response (#Delta R < 0.2) vs #eta",
                                                 resVsEtaDef.nbinsX,
                                                 resVsEtaDef.xmin,
                                                 resVsEtaDef.xmax);

  h_L1Puppi_pion_response_0p2_pt_barrel_ = ibooker.book1D("L1PUPPIPionResponse0p2VsPtBarrel",
                                                          "L1 PUPPI Pion Response (#Delta R < 0.2) vs p_{T}, Barrel",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1Puppi_pion_response_0p2_pt_endcap_ = ibooker.book1D("L1PUPPIPionResponse0p2VsPtEndcap",
                                                          "L1 PUPPI Pion Response (#Delta R < 0.2) vs p_{T}, Endcap",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1Puppi_pion_response_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIPionResponse0p2VsPtEndcapNoTk",
                     "L1 PUPPI Pion Response (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pion_response_0p2_pt_hf_ = ibooker.book1D("L1PUPPIPionResponse0p2VsPtHF",
                                                      "L1 PUPPI Pion Response (#Delta R < 0.2) vs p_{T}, HF",
                                                      resVsPtDef.nbinsX,
                                                      resVsPtDef.xmin,
                                                      resVsPtDef.xmax);

  h_L1Puppi_pion_response_0p2_eta_ = ibooker.book1D("L1PUPPIPionResponse0p2VsEta",
                                                    "L1 PUPPI Pion Response (#Delta R < 0.2) vs #eta",
                                                    resVsEtaDef.nbinsX,
                                                    resVsEtaDef.xmin,
                                                    resVsEtaDef.xmax);

  h_L1PF_pion_response_best_pt_barrel_ = ibooker.book1D("L1PFPionResponseBestVsPtBarrel",
                                                        "L1 PF Pion Response (Best) vs p_{T}, Barrel",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1PF_pion_response_best_pt_endcap_ = ibooker.book1D("L1PFPionResponseBestVsPtEndcap",
                                                        "L1 PF Pion Response (Best) vs p_{T}, Endcap",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1PF_pion_response_best_pt_ecnotk_ = ibooker.book1D("L1PFPionResponseBestVsPtEndcapNoTk",
                                                        "L1 PF Pion Response (Best) vs p_{T}, Endcap No Tk",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1PF_pion_response_best_pt_hf_ = ibooker.book1D("L1PFPionResponseBestVsPtHF",
                                                    "L1 PF Pion Response (Best) vs p_{T}, HF",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax);

  h_L1PF_pion_response_best_eta_ = ibooker.book1D("L1PFPionResponseBestVsEta",
                                                  "L1 PF Pion Response (Best) vs #eta",
                                                  resVsEtaDef.nbinsX,
                                                  resVsEtaDef.xmin,
                                                  resVsEtaDef.xmax);

  h_L1Puppi_pion_response_best_pt_barrel_ = ibooker.book1D("L1PUPPIPionResponseBestVsPtBarrel",
                                                           "L1 PUPPI Pion Response (Best) vs p_{T}, Barrel",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1Puppi_pion_response_best_pt_endcap_ = ibooker.book1D("L1PUPPIPionResponseBestVsPtEndcap",
                                                           "L1 PUPPI Pion Response (Best) vs p_{T}, Endcap",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1Puppi_pion_response_best_pt_ecnotk_ = ibooker.book1D("L1PUPPIPionResponseBestVsPtEndcapNoTk",
                                                           "L1 PUPPI Pion Response (Best) vs p_{T}, Endcap No Tk",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1Puppi_pion_response_best_pt_hf_ = ibooker.book1D("L1PUPPIPionResponseBestVsPtHF",
                                                       "L1 PUPPI Pion Response (Best) vs p_{T}, HF",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1Puppi_pion_response_best_eta_ = ibooker.book1D("L1PUPPIPionResponseBestVsEta",
                                                     "L1 PUPPI Pion Response (Best) vs #eta",
                                                     resVsEtaDef.nbinsX,
                                                     resVsEtaDef.xmin,
                                                     resVsEtaDef.xmax);

  h_L1PF_jet_response_pt_barrel_ = ibooker.book1D("L1PFJetResponseVsPtBarrel",
                                                  "L1 PF Jet Response vs p_{T}, Barrel",
                                                  resVsPtDef.nbinsX,
                                                  resVsPtDef.xmin,
                                                  resVsPtDef.xmax);

  h_L1PF_jet_response_pt_endcap_ = ibooker.book1D("L1PFJetResponseVsPtEndcap",
                                                  "L1 PF Jet Response vs p_{T}, Endcap",
                                                  resVsPtDef.nbinsX,
                                                  resVsPtDef.xmin,
                                                  resVsPtDef.xmax);

  h_L1PF_jet_response_pt_ecnotk_ = ibooker.book1D("L1PFJetResponseVsPtEndcapNoTk",
                                                  "L1 PF Jet Response vs p_{T}, Endcap No Tk",
                                                  resVsPtDef.nbinsX,
                                                  resVsPtDef.xmin,
                                                  resVsPtDef.xmax);

  h_L1PF_jet_response_pt_hf_ = ibooker.book1D(
      "L1PFJetResponseVsPtHF", "L1 PF Jet Response vs p_{T}, HF", resVsPtDef.nbinsX, resVsPtDef.xmin, resVsPtDef.xmax);

  h_L1PF_jet_response_eta_ = ibooker.book1D(
      "L1PFJetResponseVsEta", "L1 PF Jet Response vs #eta", resVsEtaDef.nbinsX, resVsEtaDef.xmin, resVsEtaDef.xmax);

  h_L1Puppi_jet_response_pt_barrel_ = ibooker.book1D("L1PUPPIJetResponseVsPtBarrel",
                                                     "L1 PUPPI Jet Response vs p_{T}, Barrel",
                                                     resVsPtDef.nbinsX,
                                                     resVsPtDef.xmin,
                                                     resVsPtDef.xmax);

  h_L1Puppi_jet_response_pt_endcap_ = ibooker.book1D("L1PUPPIJetResponseVsPtEndcap",
                                                     "L1 PUPPI Jet Response vs p_{T}, Endcap",
                                                     resVsPtDef.nbinsX,
                                                     resVsPtDef.xmin,
                                                     resVsPtDef.xmax);

  h_L1Puppi_jet_response_pt_ecnotk_ = ibooker.book1D("L1PUPPIJetResponseVsPtEndcapNoTk",
                                                     "L1 PUPPI Jet Response vs p_{T}, EndcapNoTk",
                                                     resVsPtDef.nbinsX,
                                                     resVsPtDef.xmin,
                                                     resVsPtDef.xmax);

  h_L1Puppi_jet_response_pt_hf_ = ibooker.book1D("L1PUPPIJetResponseVsPtHF",
                                                 "L1 PUPPI Jet Response vs p_{T}, HF",
                                                 resVsPtDef.nbinsX,
                                                 resVsPtDef.xmin,
                                                 resVsPtDef.xmax);

  h_L1Puppi_jet_response_eta_ = ibooker.book1D("L1PUPPIJetResponseVsEta",
                                               "L1 PUPPI Jet Response vs #eta",
                                               resVsEtaDef.nbinsX,
                                               resVsEtaDef.xmin,
                                               resVsEtaDef.xmax);

  //Resolution
  h_L1PF_electron_resolution_0p2_pt_barrel_ =
      ibooker.book1D("L1PFElectronResolution0p2VsPtBarrel",
                     "L1 PF Electron Resolution (#Delta R < 0.2) vs p_{T}, Barrel",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_electron_resolution_0p2_pt_endcap_ =
      ibooker.book1D("L1PFElectronResolution0p2VsPtEndcap",
                     "L1 PF Electron Resolution (#Delta R < 0.2) vs p_{T}, Endcap",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_electron_resolution_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PFElectronResolution0p2VsPtEndcapNoTk",
                     "L1 PF Electron Resolution (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_electron_resolution_0p2_pt_hf_ = ibooker.book1D("L1PFElectronResolution0p2VsPtHF",
                                                         "L1 PF Electron Resolution (#Delta R < 0.2) vs p_{T}, HF",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_0p2_pt_barrel_ =
      ibooker.book1D("L1PUPPIElectronResolution0p2VsPtBarrel",
                     "L1 PUPPI Electron Resolution (#Delta R < 0.2) vs p_{T}, Barrel",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_0p2_pt_endcap_ =
      ibooker.book1D("L1PUPPIElectronResolution0p2VsPtEndcap",
                     "L1 PUPPI Electron Resolution (#Delta R < 0.2) vs p_{T}, Endcap",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIElectronResolution0p2VsPtEndcapNoTk",
                     "L1 PUPPI Electron Resolution (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_0p2_pt_hf_ =
      ibooker.book1D("L1PUPPIElectronResolution0p2VsPtHF",
                     "L1 PUPPI Electron Resolution (#Delta R < 0.2) vs p_{T}, HF",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_electron_resolution_best_pt_barrel_ = ibooker.book1D("L1PFElectronResolutionBestVsPtBarrel",
                                                              "L1 PF Electron Resolution (Best) vs p_{T}, Barrel",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax);

  h_L1PF_electron_resolution_best_pt_endcap_ = ibooker.book1D("L1PFElectronResolutionBestVsPtEndcap",
                                                              "L1 PF Electron Resolution (Best) vs p_{T}, Endcap",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax);

  h_L1PF_electron_resolution_best_pt_ecnotk_ = ibooker.book1D("L1PFElectronResolutionBestVsPtEndcapNoTk",
                                                              "L1 PF Electron Resolution (Best) vs p_{T}, Endcap No Tk",
                                                              resVsPtDef.nbinsX,
                                                              resVsPtDef.xmin,
                                                              resVsPtDef.xmax);

  h_L1PF_electron_resolution_best_pt_hf_ = ibooker.book1D("L1PFElectronResolutionBestVsPtHF",
                                                          "L1 PF Electron Resolution (Best) vs p_{T}, HF",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_best_pt_barrel_ = ibooker.book1D("L1PUPPIElectronResolutionBestVsPtBarrel",
                                                                 "L1 PUPPI Electron Resolution (Best) vs p_{T}, Barrel",
                                                                 resVsPtDef.nbinsX,
                                                                 resVsPtDef.xmin,
                                                                 resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_best_pt_endcap_ = ibooker.book1D("L1PUPPIElectronResolutionBestVsPtEndcap",
                                                                 "L1 PUPPI Electron Resolution (Best) vs p_{T}, Endcap",
                                                                 resVsPtDef.nbinsX,
                                                                 resVsPtDef.xmin,
                                                                 resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_best_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIElectronResolutionBestVsPtEndcapNoTk",
                     "L1 PUPPI Electron Resolution (Best) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_electron_resolution_best_pt_hf_ = ibooker.book1D("L1PUPPIElectronResolutionBestVsPtHF",
                                                             "L1 PUPPI Electron Resolution (Best) vs p_{T}, HF",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1PF_pizero_resolution_0p2_pt_barrel_ = ibooker.book1D("L1PFPi0Resolution0p2VsPtBarrel",
                                                           "L1 PF Pi0 Resolution (#Delta R < 0.2) vs p_{T}, Barrel",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1PF_pizero_resolution_0p2_pt_endcap_ = ibooker.book1D("L1PFPi0Resolution0p2VsPtEndcap",
                                                           "L1 PF Pi0 Resolution (#Delta R < 0.2) vs p_{T}, Endcap",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1PF_pizero_resolution_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PFPi0Resolution0p2VsPtEndcapNoTk",
                     "L1 PF Pi0 Resolution (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_pizero_resolution_0p2_pt_hf_ = ibooker.book1D("L1PFPi0Resolution0p2VsPtHF",
                                                       "L1 PF Pi0 Resolution (#Delta R < 0.2) vs p_{T}, HF",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_0p2_pt_barrel_ =
      ibooker.book1D("L1PUPPIPi0Resolution0p2VsPtBarrel",
                     "L1 PUPPI Pi0 Resolution (#Delta R < 0.2) vs p_{T}, Barrel",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_0p2_pt_endcap_ =
      ibooker.book1D("L1PUPPIPi0Resolution0p2VsPtEndcap",
                     "L1 PUPPI Pi0 Resolution (#Delta R < 0.2) vs p_{T}, Endcap",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIPi0Resolution0p2VsPtEndcapNoTk",
                     "L1 PUPPI Pi0 Resolution (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_0p2_pt_hf_ = ibooker.book1D("L1PUPPIPi0Resolution0p2VsPtHF",
                                                          "L1 PUPPI Pi0 Resolution (#Delta R < 0.2) vs p_{T}, HF",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pizero_resolution_best_pt_barrel_ = ibooker.book1D("L1PFPi0ResolutionBestVsPtBarrel",
                                                            "L1 PF Pi0 Resolution (Best) vs p_{T}, Barrel",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1PF_pizero_resolution_best_pt_endcap_ = ibooker.book1D("L1PFPi0ResolutionBestVsPtEndcap",
                                                            "L1 PF Pi0 Resolution (Best) vs p_{T}, Endcap",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1PF_pizero_resolution_best_pt_ecnotk_ = ibooker.book1D("L1PFPi0ResolutionBestVsPtEndcapNoTk",
                                                            "L1 PF Pi0 Resolution (Best) vs p_{T}, Endcap No Tk",
                                                            resVsPtDef.nbinsX,
                                                            resVsPtDef.xmin,
                                                            resVsPtDef.xmax);

  h_L1PF_pizero_resolution_best_pt_hf_ = ibooker.book1D("L1PFPi0ResolutionBestVsPtHF",
                                                        "L1 PF Pi0 Resolution (Best) vs p_{T}, HF",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_best_pt_barrel_ = ibooker.book1D("L1PUPPIPi0ResolutionBestVsPtBarrel",
                                                               "L1 PUPPI Pi0 Resolution (Best) vs p_{T}, Barrel",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_best_pt_endcap_ = ibooker.book1D("L1PUPPIPi0ResolutionBestVsPtEndcap",
                                                               "L1 PUPPI Pi0 Resolution (Best) vs p_{T}, Endcap",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_best_pt_ecnotk_ = ibooker.book1D("L1PUPPIPi0ResolutionBestVsPtEndcapNoTk",
                                                               "L1 PUPPI Pi0 Resolution (Best) vs p_{T}, Endcap No Tk",
                                                               resVsPtDef.nbinsX,
                                                               resVsPtDef.xmin,
                                                               resVsPtDef.xmax);

  h_L1Puppi_pizero_resolution_best_pt_hf_ = ibooker.book1D("L1PUPPIPi0ResolutionBestVsPtHF",
                                                           "L1 PUPPI Pi0 Resolution (Best) vs p_{T}, HF",
                                                           resVsPtDef.nbinsX,
                                                           resVsPtDef.xmin,
                                                           resVsPtDef.xmax);

  h_L1PF_pion_resolution_0p2_pt_barrel_ = ibooker.book1D("L1PFPionResolution0p2VsPtBarrel",
                                                         "L1 PF Pion Resolution (#Delta R < 0.2) vs p_{T}, Barrel",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1PF_pion_resolution_0p2_pt_endcap_ = ibooker.book1D("L1PFPionResolution0p2VsPtEndcap",
                                                         "L1 PF Pion Resolution (#Delta R < 0.2) vs p_{T}, Endcap",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1PF_pion_resolution_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PFPionResolution0p2VsPtEndcapNoTk",
                     "L1 PF Pion Resolution (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1PF_pion_resolution_0p2_pt_hf_ = ibooker.book1D("L1PFPionResolution0p2VsPtHF",
                                                     "L1 PF Pion Resolution (#Delta R < 0.2) vs p_{T}, HF",
                                                     resVsPtDef.nbinsX,
                                                     resVsPtDef.xmin,
                                                     resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_0p2_pt_barrel_ =
      ibooker.book1D("L1PUPPIPionResolution0p2VsPtBarrel",
                     "L1 PUPPI Pion Resolution (#Delta R < 0.2) vs p_{T}, Barrel",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_0p2_pt_endcap_ =
      ibooker.book1D("L1PUPPIPionResolution0p2VsPtEndcap",
                     "L1 PUPPI Pion Resolution (#Delta R < 0.2) vs p_{T}, Endcap",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_0p2_pt_ecnotk_ =
      ibooker.book1D("L1PUPPIPionResolution0p2VsPtEndcapNoTk",
                     "L1 PUPPI Pion Resolution (#Delta R < 0.2) vs p_{T}, Endcap No Tk",
                     resVsPtDef.nbinsX,
                     resVsPtDef.xmin,
                     resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_0p2_pt_hf_ = ibooker.book1D("L1PUPPIPionResolution0p2VsPtHF",
                                                        "L1 PUPPI Pion Resolution (#Delta R < 0.2) vs p_{T}, HF",
                                                        resVsPtDef.nbinsX,
                                                        resVsPtDef.xmin,
                                                        resVsPtDef.xmax);

  h_L1PF_pion_resolution_best_pt_barrel_ = ibooker.book1D("L1PFPionResolutionBestVsPtBarrel",
                                                          "L1 PF Pion Resolution (Best) vs p_{T}, Barrel",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pion_resolution_best_pt_endcap_ = ibooker.book1D("L1PFPionResolutionBestVsPtEndcap",
                                                          "L1 PF Pion Resolution (Best) vs p_{T}, Endcap",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pion_resolution_best_pt_ecnotk_ = ibooker.book1D("L1PFPionResolutionBestVsPtEndcapNoTk",
                                                          "L1 PF Pion Resolution (Best) vs p_{T}, Endcap No Tk",
                                                          resVsPtDef.nbinsX,
                                                          resVsPtDef.xmin,
                                                          resVsPtDef.xmax);

  h_L1PF_pion_resolution_best_pt_hf_ = ibooker.book1D("L1PFPionResolutionBestVsPtHF",
                                                      "L1 PF Pion Resolution (Best) vs p_{T}, HF",
                                                      resVsPtDef.nbinsX,
                                                      resVsPtDef.xmin,
                                                      resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_best_pt_barrel_ = ibooker.book1D("L1PUPPIPionResolutionBestVsPtBarrel",
                                                             "L1 PUPPI Pion Resolution (Best) vs p_{T}, Barrel",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_best_pt_endcap_ = ibooker.book1D("L1PUPPIPionResolutionBestVsPtEndcap",
                                                             "L1 PUPPI Pion Resolution (Best) vs p_{T}, Endcap",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_best_pt_ecnotk_ = ibooker.book1D("L1PUPPIPionResolutionBestVsPtEndcapNoTk",
                                                             "L1 PUPPI Pion Resolution (Best) vs p_{T}, Endcap No Tk",
                                                             resVsPtDef.nbinsX,
                                                             resVsPtDef.xmin,
                                                             resVsPtDef.xmax);

  h_L1Puppi_pion_resolution_best_pt_hf_ = ibooker.book1D("L1PUPPIPionResolutionBestVsPtHF",
                                                         "L1 PUPPI Pion Resolution (Best) vs p_{T}, HF",
                                                         resVsPtDef.nbinsX,
                                                         resVsPtDef.xmin,
                                                         resVsPtDef.xmax);

  h_L1PF_jet_resolution_pt_barrel_ = ibooker.book1D("L1PFJetResolutionVsPtBarrel",
                                                    "L1 PF Jet Resolution vs p_{T}, Barrel",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax);

  h_L1PF_jet_resolution_pt_endcap_ = ibooker.book1D("L1PFJetResolutionVsPtEndcap",
                                                    "L1 PF Jet Resolution vs p_{T}, Endcap",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax);

  h_L1PF_jet_resolution_pt_ecnotk_ = ibooker.book1D("L1PFJetResolutionVsPtEndcapNoTk",
                                                    "L1 PF Jet Resolution vs p_{T}, Endcap No Tk",
                                                    resVsPtDef.nbinsX,
                                                    resVsPtDef.xmin,
                                                    resVsPtDef.xmax);

  h_L1PF_jet_resolution_pt_hf_ = ibooker.book1D("L1PFJetResolutionVsPtHF",
                                                "L1 PF Jet Resolution vs p_{T}, HF",
                                                resVsPtDef.nbinsX,
                                                resVsPtDef.xmin,
                                                resVsPtDef.xmax);

  h_L1Puppi_jet_resolution_pt_barrel_ = ibooker.book1D("L1PUPPIJetResolutionVsPtBarrel",
                                                       "L1 PUPPI Jet Resolution vs p_{T}, Barrel",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1Puppi_jet_resolution_pt_endcap_ = ibooker.book1D("L1PUPPIJetResolutionVsPtEndcap",
                                                       "L1 PUPPI Jet Resolution vs p_{T}, Endcap",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1Puppi_jet_resolution_pt_ecnotk_ = ibooker.book1D("L1PUPPIJetResolutionVsPtEndcapNoTk",
                                                       "L1 PUPPI Jet Resolution vs p_{T}, EndcapNoTk",
                                                       resVsPtDef.nbinsX,
                                                       resVsPtDef.xmin,
                                                       resVsPtDef.xmax);

  h_L1Puppi_jet_resolution_pt_hf_ = ibooker.book1D("L1PUPPIJetResolutionVsPtHF",
                                                   "L1 PUPPI Jet Resolution vs p_{T}, HF",
                                                   resVsPtDef.nbinsX,
                                                   resVsPtDef.xmin,
                                                   resVsPtDef.xmax);

  ibooker.cd();

  return;
}

//
// -------------------------------------- endRun --------------------------------------------
//
void L1TPhase2CorrelatorOffline::dqmEndRun(const edm::Run& run, const edm::EventSetup& iSetup) {
  computeResponseResolution();
}

void L1TPhase2CorrelatorOffline::computeResponseResolution() {
  std::vector<MonitorElement*> monElementstoComputeIn = {h_L1PF_electron_ptratio_0p2_vs_pt_barrel_,
                                                         h_L1PF_electron_ptratio_0p2_vs_pt_endcap_,
                                                         h_L1PF_electron_ptratio_0p2_vs_pt_ecnotk_,
                                                         h_L1PF_electron_ptratio_0p2_vs_pt_hf_,
                                                         h_L1PF_electron_ptratio_0p2_vs_eta_,
                                                         h_L1Puppi_electron_ptratio_0p2_vs_pt_barrel_,
                                                         h_L1Puppi_electron_ptratio_0p2_vs_pt_endcap_,
                                                         h_L1Puppi_electron_ptratio_0p2_vs_pt_ecnotk_,
                                                         h_L1Puppi_electron_ptratio_0p2_vs_pt_hf_,
                                                         h_L1Puppi_electron_ptratio_0p2_vs_eta_,
                                                         h_L1PF_electron_ptratio_best_vs_pt_barrel_,
                                                         h_L1PF_electron_ptratio_best_vs_pt_endcap_,
                                                         h_L1PF_electron_ptratio_best_vs_pt_ecnotk_,
                                                         h_L1PF_electron_ptratio_best_vs_pt_hf_,
                                                         h_L1PF_electron_ptratio_best_vs_eta_,
                                                         h_L1Puppi_electron_ptratio_best_vs_pt_barrel_,
                                                         h_L1Puppi_electron_ptratio_best_vs_pt_endcap_,
                                                         h_L1Puppi_electron_ptratio_best_vs_pt_ecnotk_,
                                                         h_L1Puppi_electron_ptratio_best_vs_pt_hf_,
                                                         h_L1Puppi_electron_ptratio_best_vs_eta_,
                                                         h_L1PF_pizero_ptratio_0p2_vs_pt_barrel_,
                                                         h_L1PF_pizero_ptratio_0p2_vs_pt_endcap_,
                                                         h_L1PF_pizero_ptratio_0p2_vs_pt_ecnotk_,
                                                         h_L1PF_pizero_ptratio_0p2_vs_pt_hf_,
                                                         h_L1PF_pizero_ptratio_0p2_vs_eta_,
                                                         h_L1Puppi_pizero_ptratio_0p2_vs_pt_barrel_,
                                                         h_L1Puppi_pizero_ptratio_0p2_vs_pt_endcap_,
                                                         h_L1Puppi_pizero_ptratio_0p2_vs_pt_ecnotk_,
                                                         h_L1Puppi_pizero_ptratio_0p2_vs_pt_hf_,
                                                         h_L1Puppi_pizero_ptratio_0p2_vs_eta_,
                                                         h_L1PF_pizero_ptratio_best_vs_pt_barrel_,
                                                         h_L1PF_pizero_ptratio_best_vs_pt_endcap_,
                                                         h_L1PF_pizero_ptratio_best_vs_pt_ecnotk_,
                                                         h_L1PF_pizero_ptratio_best_vs_pt_hf_,
                                                         h_L1PF_pizero_ptratio_best_vs_eta_,
                                                         h_L1Puppi_pizero_ptratio_best_vs_pt_barrel_,
                                                         h_L1Puppi_pizero_ptratio_best_vs_pt_endcap_,
                                                         h_L1Puppi_pizero_ptratio_best_vs_pt_ecnotk_,
                                                         h_L1Puppi_pizero_ptratio_best_vs_pt_hf_,
                                                         h_L1Puppi_pizero_ptratio_best_vs_eta_,
                                                         h_L1PF_pion_ptratio_0p2_vs_pt_barrel_,
                                                         h_L1PF_pion_ptratio_0p2_vs_pt_endcap_,
                                                         h_L1PF_pion_ptratio_0p2_vs_pt_ecnotk_,
                                                         h_L1PF_pion_ptratio_0p2_vs_pt_hf_,
                                                         h_L1PF_pion_ptratio_0p2_vs_eta_,
                                                         h_L1Puppi_pion_ptratio_0p2_vs_pt_barrel_,
                                                         h_L1Puppi_pion_ptratio_0p2_vs_pt_endcap_,
                                                         h_L1Puppi_pion_ptratio_0p2_vs_pt_ecnotk_,
                                                         h_L1Puppi_pion_ptratio_0p2_vs_pt_hf_,
                                                         h_L1Puppi_pion_ptratio_0p2_vs_eta_,
                                                         h_L1PF_pion_ptratio_best_vs_pt_barrel_,
                                                         h_L1PF_pion_ptratio_best_vs_pt_endcap_,
                                                         h_L1PF_pion_ptratio_best_vs_pt_ecnotk_,
                                                         h_L1PF_pion_ptratio_best_vs_pt_hf_,
                                                         h_L1PF_pion_ptratio_best_vs_eta_,
                                                         h_L1Puppi_pion_ptratio_best_vs_pt_barrel_,
                                                         h_L1Puppi_pion_ptratio_best_vs_pt_endcap_,
                                                         h_L1Puppi_pion_ptratio_best_vs_pt_ecnotk_,
                                                         h_L1Puppi_pion_ptratio_best_vs_pt_hf_,
                                                         h_L1Puppi_pion_ptratio_best_vs_eta_,
                                                         h_L1PF_jet_ptratio_vs_pt_barrel_,
                                                         h_L1PF_jet_ptratio_vs_pt_endcap_,
                                                         h_L1PF_jet_ptratio_vs_pt_ecnotk_,
                                                         h_L1PF_jet_ptratio_vs_pt_hf_,
                                                         h_L1PF_jet_ptratio_vs_eta_,
                                                         h_L1Puppi_jet_ptratio_vs_pt_barrel_,
                                                         h_L1Puppi_jet_ptratio_vs_pt_endcap_,
                                                         h_L1Puppi_jet_ptratio_vs_pt_ecnotk_,
                                                         h_L1Puppi_jet_ptratio_vs_pt_hf_,
                                                         h_L1Puppi_jet_ptratio_vs_eta_};
  std::vector<MonitorElement*> monElementstoComputeResp = {h_L1PF_electron_response_0p2_pt_barrel_,
                                                           h_L1PF_electron_response_0p2_pt_endcap_,
                                                           h_L1PF_electron_response_0p2_pt_ecnotk_,
                                                           h_L1PF_electron_response_0p2_pt_hf_,
                                                           h_L1PF_electron_response_0p2_eta_,
                                                           h_L1Puppi_electron_response_0p2_pt_barrel_,
                                                           h_L1Puppi_electron_response_0p2_pt_endcap_,
                                                           h_L1Puppi_electron_response_0p2_pt_ecnotk_,
                                                           h_L1Puppi_electron_response_0p2_pt_hf_,
                                                           h_L1Puppi_electron_response_0p2_eta_,
                                                           h_L1PF_electron_response_best_pt_barrel_,
                                                           h_L1PF_electron_response_best_pt_endcap_,
                                                           h_L1PF_electron_response_best_pt_ecnotk_,
                                                           h_L1PF_electron_response_best_pt_hf_,
                                                           h_L1PF_electron_response_best_eta_,
                                                           h_L1Puppi_electron_response_best_pt_barrel_,
                                                           h_L1Puppi_electron_response_best_pt_endcap_,
                                                           h_L1Puppi_electron_response_best_pt_ecnotk_,
                                                           h_L1Puppi_electron_response_best_pt_hf_,
                                                           h_L1Puppi_electron_response_best_eta_,
                                                           h_L1PF_pizero_response_0p2_pt_barrel_,
                                                           h_L1PF_pizero_response_0p2_pt_endcap_,
                                                           h_L1PF_pizero_response_0p2_pt_ecnotk_,
                                                           h_L1PF_pizero_response_0p2_pt_hf_,
                                                           h_L1PF_pizero_response_0p2_eta_,
                                                           h_L1Puppi_pizero_response_0p2_pt_barrel_,
                                                           h_L1Puppi_pizero_response_0p2_pt_endcap_,
                                                           h_L1Puppi_pizero_response_0p2_pt_ecnotk_,
                                                           h_L1Puppi_pizero_response_0p2_pt_hf_,
                                                           h_L1Puppi_pizero_response_0p2_eta_,
                                                           h_L1PF_pizero_response_best_pt_barrel_,
                                                           h_L1PF_pizero_response_best_pt_endcap_,
                                                           h_L1PF_pizero_response_best_pt_ecnotk_,
                                                           h_L1PF_pizero_response_best_pt_hf_,
                                                           h_L1PF_pizero_response_best_eta_,
                                                           h_L1Puppi_pizero_response_best_pt_barrel_,
                                                           h_L1Puppi_pizero_response_best_pt_endcap_,
                                                           h_L1Puppi_pizero_response_best_pt_ecnotk_,
                                                           h_L1Puppi_pizero_response_best_pt_hf_,
                                                           h_L1Puppi_pizero_response_best_eta_,
                                                           h_L1PF_pion_response_0p2_pt_barrel_,
                                                           h_L1PF_pion_response_0p2_pt_endcap_,
                                                           h_L1PF_pion_response_0p2_pt_ecnotk_,
                                                           h_L1PF_pion_response_0p2_pt_hf_,
                                                           h_L1PF_pion_response_0p2_eta_,
                                                           h_L1Puppi_pion_response_0p2_pt_barrel_,
                                                           h_L1Puppi_pion_response_0p2_pt_endcap_,
                                                           h_L1Puppi_pion_response_0p2_pt_ecnotk_,
                                                           h_L1Puppi_pion_response_0p2_pt_hf_,
                                                           h_L1Puppi_pion_response_0p2_eta_,
                                                           h_L1PF_pion_response_best_pt_barrel_,
                                                           h_L1PF_pion_response_best_pt_endcap_,
                                                           h_L1PF_pion_response_best_pt_ecnotk_,
                                                           h_L1PF_pion_response_best_pt_hf_,
                                                           h_L1PF_pion_response_best_eta_,
                                                           h_L1Puppi_pion_response_best_pt_barrel_,
                                                           h_L1Puppi_pion_response_best_pt_endcap_,
                                                           h_L1Puppi_pion_response_best_pt_ecnotk_,
                                                           h_L1Puppi_pion_response_best_pt_hf_,
                                                           h_L1Puppi_pion_response_best_eta_,
                                                           h_L1PF_jet_response_pt_barrel_,
                                                           h_L1PF_jet_response_pt_endcap_,
                                                           h_L1PF_jet_response_pt_ecnotk_,
                                                           h_L1PF_jet_response_pt_hf_,
                                                           h_L1PF_jet_response_eta_,
                                                           h_L1Puppi_jet_response_pt_barrel_,
                                                           h_L1Puppi_jet_response_pt_endcap_,
                                                           h_L1Puppi_jet_response_pt_ecnotk_,
                                                           h_L1Puppi_jet_response_pt_hf_,
                                                           h_L1Puppi_jet_response_eta_};
  std::vector<MonitorElement*> monElementstoComputeResol = {
      h_L1PF_electron_resolution_0p2_pt_barrel_,
      h_L1PF_electron_resolution_0p2_pt_endcap_,
      h_L1PF_electron_resolution_0p2_pt_ecnotk_,
      h_L1PF_electron_resolution_0p2_pt_hf_,
      nullptr, 
      h_L1Puppi_electron_resolution_0p2_pt_barrel_,
      h_L1Puppi_electron_resolution_0p2_pt_endcap_,
      h_L1Puppi_electron_resolution_0p2_pt_ecnotk_,
      h_L1Puppi_electron_resolution_0p2_pt_hf_,
      nullptr, 
      h_L1PF_electron_resolution_best_pt_barrel_,
      h_L1PF_electron_resolution_best_pt_endcap_,
      h_L1PF_electron_resolution_best_pt_ecnotk_,
      h_L1PF_electron_resolution_best_pt_hf_,
      nullptr, 
      h_L1Puppi_electron_resolution_best_pt_barrel_,
      h_L1Puppi_electron_resolution_best_pt_endcap_,
      h_L1Puppi_electron_resolution_best_pt_ecnotk_,
      h_L1Puppi_electron_resolution_best_pt_hf_,
      nullptr, 
      h_L1PF_pizero_resolution_0p2_pt_barrel_,
      h_L1PF_pizero_resolution_0p2_pt_endcap_,
      h_L1PF_pizero_resolution_0p2_pt_ecnotk_,
      h_L1PF_pizero_resolution_0p2_pt_hf_,
      nullptr, 
      h_L1Puppi_pizero_resolution_0p2_pt_barrel_,
      h_L1Puppi_pizero_resolution_0p2_pt_endcap_,
      h_L1Puppi_pizero_resolution_0p2_pt_ecnotk_,
      h_L1Puppi_pizero_resolution_0p2_pt_hf_,
      nullptr, 
      h_L1PF_pizero_resolution_best_pt_barrel_,
      h_L1PF_pizero_resolution_best_pt_endcap_,
      h_L1PF_pizero_resolution_best_pt_ecnotk_,
      h_L1PF_pizero_resolution_best_pt_hf_,
      nullptr, 
      h_L1Puppi_pizero_resolution_best_pt_barrel_,
      h_L1Puppi_pizero_resolution_best_pt_endcap_,
      h_L1Puppi_pizero_resolution_best_pt_ecnotk_,
      h_L1Puppi_pizero_resolution_best_pt_hf_,
      nullptr, 
      h_L1PF_pion_resolution_0p2_pt_barrel_,
      h_L1PF_pion_resolution_0p2_pt_endcap_,
      h_L1PF_pion_resolution_0p2_pt_ecnotk_,
      h_L1PF_pion_resolution_0p2_pt_hf_,
      nullptr, 
      h_L1Puppi_pion_resolution_0p2_pt_barrel_,
      h_L1Puppi_pion_resolution_0p2_pt_endcap_,
      h_L1Puppi_pion_resolution_0p2_pt_ecnotk_,
      h_L1Puppi_pion_resolution_0p2_pt_hf_,
      nullptr, 
      h_L1PF_pion_resolution_best_pt_barrel_,
      h_L1PF_pion_resolution_best_pt_endcap_,
      h_L1PF_pion_resolution_best_pt_ecnotk_,
      h_L1PF_pion_resolution_best_pt_hf_,
      nullptr, 
      h_L1Puppi_pion_resolution_best_pt_barrel_,
      h_L1Puppi_pion_resolution_best_pt_endcap_,
      h_L1Puppi_pion_resolution_best_pt_ecnotk_,
      h_L1Puppi_pion_resolution_best_pt_hf_,
      nullptr, 
      h_L1PF_jet_resolution_pt_barrel_,
      h_L1PF_jet_resolution_pt_endcap_,
      h_L1PF_jet_resolution_pt_ecnotk_,
      h_L1PF_jet_resolution_pt_hf_,
      nullptr, 
      h_L1Puppi_jet_resolution_pt_barrel_,
      h_L1Puppi_jet_resolution_pt_endcap_,
      h_L1Puppi_jet_resolution_pt_ecnotk_,
      h_L1Puppi_jet_resolution_pt_hf_,
      nullptr 
  };

  for (unsigned int i = 0; i < monElementstoComputeIn.size(); i++) {
    if (monElementstoComputeIn[i] != nullptr && monElementstoComputeResp[i] != nullptr &&
        monElementstoComputeResol[i] != nullptr) {
      medianResponseCorrResolution(
          monElementstoComputeIn[i], monElementstoComputeResp[i], monElementstoComputeResol[i]);
    } else if (monElementstoComputeIn[i] != nullptr && monElementstoComputeResp[i] != nullptr) {
      medianResponse(monElementstoComputeIn[i], monElementstoComputeResp[i]);
    }
  }
}

void L1TPhase2CorrelatorOffline::medianResponse(MonitorElement* in2D, MonitorElement* response) {
  auto hbase = in2D->getTH2F();
  auto hresp = response->getTH1F();
  if (hbase != nullptr && hresp != nullptr) {
    if (hbase->GetNbinsX() == hresp->GetNbinsX() && hbase->GetNbinsX()) {
      auto med = hbase->QuantilesX(0.5, "_qx");
      for (int ib = 0; ib < hbase->GetNbinsX() + 1; ib++) {
        hresp->SetBinContent(ib, med->GetBinContent(ib));
      }
    }
  }
}

void L1TPhase2CorrelatorOffline::medianResponseCorrResolution(MonitorElement* in2D,
                                                              MonitorElement* response,
                                                              MonitorElement* resolution) {
  auto hbase = in2D->getTH2F();
  auto hresp = response->getTH1F();
  auto hresol = resolution->getTH1F();
  if (hbase != nullptr && hresp != nullptr && hresol != nullptr) {
    if (hbase->GetNbinsX() == hresp->GetNbinsX() && hbase->GetNbinsX() == hresol->GetNbinsX()) {
      auto med = hbase->QuantilesX(0.5, "_qx");
      TGraph* ptrecgen = new TGraph(hbase->GetNbinsX());
      for (int ib = 1; ib < hbase->GetNbinsX() + 1; ib++) {
        float corr = med->GetBinContent(ib);
        float xval = hbase->GetXaxis()->GetBinCenter(ib);
        ptrecgen->SetPoint(ib - 1, xval * corr, xval);
        hresp->SetBinContent(ib, corr);
      }
      delete med;
      ptrecgen->Sort();
      TH2F* ch = new TH2F(*hbase);
      ch->Reset("ICE");
      for (int ibx = 1; ibx < ch->GetNbinsX() + 1; ibx++) {
        float xval = hbase->GetXaxis()->GetBinCenter(ibx);
        for (int iby = 1; iby < ch->GetNbinsY() + 1; iby++) {
          float yval = hbase->GetYaxis()->GetBinCenter(iby);
          float newyval = ptrecgen->Eval(yval * xval) / xval;
          int ycb = ch->FindBin(xval, newyval);
          ch->SetBinContent(ycb, ch->GetBinContent(ycb) + hbase->GetBinContent(ibx, iby));
        }
      }
      delete ptrecgen;
      auto qc = ch->QuantilesX(0.5, "_qc");
      auto qhi = ch->QuantilesX(0.84, "_qhi");
      auto qlo = ch->QuantilesX(0.16, "_qlo");
      delete ch;
      for (int ibx = 1; ibx < hbase->GetNbinsX() + 1; ibx++) {
        hresol->SetBinContent(
            ibx, qc->GetBinContent(ibx) > 0.2 ? (qhi->GetBinContent(ibx) - qlo->GetBinContent(ibx)) / 2. : 0.);
      }
      delete qc;
      delete qhi;
      delete qlo;
    }
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(L1TPhase2CorrelatorOffline);
