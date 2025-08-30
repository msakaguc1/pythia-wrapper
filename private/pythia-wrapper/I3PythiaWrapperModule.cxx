// I3PythiaWrapperModule.cxx
// This file creates an I3MCTree in various formats by correctly ordering outputs from Pythia 
// Any functions/calls can be deleted with the exception of the first three passes. Keys created from this file are:
////HadronizationTree: Saves all the immediate children as a product of hadronization from the chain of quarks and the muon/numu from CC or NC DIS
////I3MCTree: The tree of all "stable particles" attached to the primary neutrino in accordance to Pythia's definition and default lifetime
////HadronTypedMCTree: Saves all hadrons as type "Hadron" instead of their individual
////ParamCascade: Hadron Blob stored in as a I3Particle created by MakeParamCascade
////ParamTree: The hadron blob stored with the remaining particles, created by MakeParamMCTree
////IntermediateMCTree: an I3MCTree with all intermediate particles (excluding the gluons, including the quarks)
////NCRecoilDarkTree: Darkened the reociling neutrino in Neutral Current DIS
#include "pythia-wrapper/I3PythiaWrapperModule.hh"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/LHAPDF6.h"

#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3EventHeader.h>
#include <dataclasses/I3Constants.h>
#include <boost/make_shared.hpp>
#include <vector>
#include <cmath>
#include <map>
#include <fstream>
#include <queue>


namespace {
    // --- 1st pass: remove all intermediate gluons ---
    void PruneIntermediateGluons(I3MCTree& tree) {
        using namespace I3MCTreeUtils;
        std::vector<I3ParticleID> gluonIDs;
        for (auto it = tree.begin(); it != tree.end(); ++it)
            if (std::abs(it->GetPdgEncoding()) == 21)
                gluonIDs.push_back(it->GetID());
        for (auto const& gid : gluonIDs) {
            tree.flatten(gid);
            tree.erase(gid);
        }
    }
    // --- 2nd pass: isolate the proton branch and keep only the e_neutrino branch ---
    void RemoveProtonBranch(I3MCTree& tree) {
        using namespace I3MCTreeUtils;
        auto heads = tree.get_heads();
        for (auto const& root : heads) {
            if (root.GetPdgEncoding() == 2212) {
                log_info("Removing proton branch: PDG=%d  E=%.2f",
                     root.GetPdgEncoding(), root.GetEnergy());
                tree.erase(root.GetID());
                break;  // assume just one proton root
            }
        }
    }
    // --- 3rd pass: Delete all intermediate particles and output a MCTree with the primary neutrino and it's children
    void KeepOnlyChildrenParticles(I3MCTree& tree) {
        using namespace I3MCTreeUtils;
        std::set<I3ParticleID> headIDs;
        for (auto const& h : tree.get_heads())
            headIDs.insert(h.GetID());
        std::vector<I3ParticleID> toRemove;
        for (auto it = tree.begin(); it != tree.end(); ++it) {
            I3ParticleID pid = it->GetID();
            if (headIDs.count(pid)) continue;
            if (!tree.children(pid).empty())
                toRemove.push_back(pid);
        }
        for (auto const& pid : toRemove) {
            tree.flatten(pid);
            tree.erase(pid);
        }
    }
    void MakeHadronizationTree(I3FramePtr frame) {
    I3MCTreeConstPtr orig = frame->Get<I3MCTreeConstPtr>("I3MCTree");
    if (!orig) return;
    
    I3MCTreePtr hadronTree(new I3MCTree(*orig)); // deep copy
    
    // Find primary neutrino and direct hadronization products
    std::set<I3ParticleID> keepIDs;
    std::set<int> quark_pdg_codes = {1, 2, 3, 4, 5, 6};
    std::set<int> neutrino_pdg_codes = {12, 14, 16};
    std::set<int> lepton_pdg_codes = {11, 13, 15};
    
    using namespace I3MCTreeUtils;
    
    for (auto it = hadronTree->begin(); it != hadronTree->end(); ++it) {
        int abs_pdg = std::abs(it->GetPdgEncoding());
        I3ParticleID pid = it->GetID();
        
        // Finds primary particle
        auto primaries = GetPrimaries(*hadronTree);
        bool is_primary = false;
        for (const auto& primary : primaries) {
            if (primary.GetID() == pid && neutrino_pdg_codes.count(abs_pdg)) {
                is_primary = true;
                break;
            }
        }
        
        if (is_primary) {
            keepIDs.insert(pid);
            continue;
        }
        
        // Finds direct products of quarks (hadronization products)
        bool has_quark_parent = false;
        if (HasParent(*hadronTree, pid)) {
            try {
                I3Particle parent = GetParent(*hadronTree, pid);
                int parent_abs_pdg = std::abs(parent.GetPdgEncoding());
                if (quark_pdg_codes.count(parent_abs_pdg)) {
                    has_quark_parent = true;
                }
            } catch (...) {
                // Parent not found or other issue, skip
            }
        }
        
        if (has_quark_parent) {
            keepIDs.insert(pid);
            continue;
        }
        
        // Finds terminal leptons (muons, electrons, taus) from neutrino interaction
        if (lepton_pdg_codes.count(abs_pdg)) {
            // Check if this particle has children
            auto children = GetDaughters(*hadronTree, pid);
            bool has_children = !children.empty();
            
            // Only consider terminal particles (no children)
            if (!has_children) {
                keepIDs.insert(pid);
            }
        }
        
        // Finds terminal outgoing neutrinos from NC interactions
        if (neutrino_pdg_codes.count(abs_pdg) && !is_primary) {
            // Check if this particle has children
            auto children = GetDaughters(*hadronTree, pid);
            bool has_children = !children.empty();
            
            // Only consider terminal neutrinos with neutrino parents
            if (!has_children) {
                bool has_neutrino_parent = false;
                if (HasParent(*hadronTree, pid)) {
                    try {
                        I3Particle parent = GetParent(*hadronTree, pid);
                        int parent_abs_pdg = std::abs(parent.GetPdgEncoding());
                        if (neutrino_pdg_codes.count(parent_abs_pdg)) {
                            has_neutrino_parent = true;
                        }
                    } catch (...) {
                        // Parent not found, skip
                    }
                }
                
                if (has_neutrino_parent) {
                    keepIDs.insert(pid);
                }
            }
        }
    }
    
    // Remove any quarks from keepIDs
    std::set<I3ParticleID> finalKeepIDs;
    for (const auto& keepID : keepIDs) {
        try {
            I3Particle particle = GetParticle(*hadronTree, keepID);
            int abs_pdg = std::abs(particle.GetPdgEncoding());
            if (!quark_pdg_codes.count(abs_pdg) || lepton_pdg_codes.count(abs_pdg)) { // Only keep non-quarks
                finalKeepIDs.insert(keepID);
            }
        } catch (...) {
            continue;
        }
    }
    
    // Build clean tree with primary neutrino as root
    I3MCTreePtr cleanTree(new I3MCTree());
    
    // Find primary neutrino to use as root
    I3ParticleID primaryNuID;
    bool foundPrimaryNu = false;
    
    for (const auto& keepID : finalKeepIDs) {
        try {
            I3Particle particle = GetParticle(*hadronTree, keepID);
            int abs_pdg = std::abs(particle.GetPdgEncoding());
            if (neutrino_pdg_codes.count(abs_pdg)) {
                auto primaries = GetPrimaries(*hadronTree);
                for (const auto& primary : primaries) {
                    if (primary.GetID() == keepID) {
                        primaryNuID = keepID;
                        foundPrimaryNu = true;
                        AddPrimary(*cleanTree, particle);
                        break;
                    }
                }
                if (foundPrimaryNu) break;
            }
        } catch (...) {
            continue;
        }
    }
    
    if (!foundPrimaryNu) {
        log_warn("No primary neutrino found in HadronizationTree");
        return;
    }
    
    // Add all other kept particles as children of the primary neutrino
    for (const auto& keepID : finalKeepIDs) {
        if (keepID == primaryNuID) continue; // Skip primary neutrino itself
        
        try {
            I3Particle particle = GetParticle(*hadronTree, keepID);
            AppendChild(*cleanTree, primaryNuID, particle);
        } catch (...) {
            log_warn("Could not add particle to HadronizationTree");
        }
    }
    
    log_info("HadronizationTree: kept %zu particles", finalKeepIDs.size());
    
    frame->Put("HadronizationTree", cleanTree);
}
// --- Setting up the "hadronic blob" to replicate NuGen outputs with Pythia for CMC 
    // --- Sums all secondary hadrons into one "ParamCascade" I3Particle ---
    void MakeParamCascade(I3FramePtr frame) {
        I3MCTreeConstPtr tree = frame->Get<I3MCTreeConstPtr>("I3MCTree");
        if (!tree) return;
        std::vector<const I3Particle*> parts;
        for (const auto &p : *tree) {
            int ap = std::abs(p.GetPdgEncoding());
            if (ap!=11 && ap!=12 && ap!=13 && ap!=14 && ap!=15 && ap!=16 && ap!=22) // everything except e,μ,τ,ν and γ
                parts.push_back(&p);
        }
        // initiate I3Particle Properties
        double E_sum = 0.0;
        I3Position mom_sum(0,0,0), pos_wsum(0,0,0);
        std::vector<double> times;
        times.reserve(parts.size());

        for (auto *p : parts) {
            double e = p->GetEnergy();
            E_sum += e;
            I3Direction dir = p->GetDir();
            mom_sum  += I3Position(dir.GetX()*e, dir.GetY()*e, dir.GetZ()*e);
            
            I3Position pos = p->GetPos();
            pos_wsum += I3Position(pos.GetX()*e, pos.GetY()*e, pos.GetZ()*e);
            times.push_back(p->GetTime());
        }
        I3Particle cascade;
        cascade.SetShape(I3Particle::Cascade);
        cascade.SetEnergy(E_sum);
        

        // normalize & set direction
        if (mom_sum.Magnitude() > 0.0) {
            I3Direction direction(mom_sum.GetX(), mom_sum.GetY(), mom_sum.GetZ());
            cascade.SetDir(I3Direction(mom_sum.GetX(), mom_sum.GetY(), mom_sum.GetZ()));
        }
        // set time = earliest hadron time (or 0)
        double tmin = times.empty() ? 0.0 : *std::min_element(times.begin(), times.end());
        cascade.SetTime(tmin);

        // set position = energy‐weighted centroid
        if (E_sum > 0.0) {
            cascade.SetPos(I3Position(pos_wsum.GetX()/E_sum, pos_wsum.GetY()/E_sum, pos_wsum.GetZ()/E_sum));
        }
        cascade.SetLocationType(I3Particle::InIce);
        cascade.SetType(I3Particle::Hadrons);
        //  publish into the frame
        frame->Put("ParamCascade", boost::make_shared<I3Particle>(cascade));
    }
    void MakeParamMCTree(I3FramePtr frame) {
        I3MCTreeConstPtr orig = frame->Get<I3MCTreeConstPtr>("I3MCTree");
        I3ParticleConstPtr cascade = frame->Get<I3ParticleConstPtr>("ParamCascade");
        if (!orig || !cascade) return;

        std::vector<const I3Particle*> nus;
        for (auto &p : *orig) {
            int ap = std::abs(p.GetPdgEncoding());
            if (ap==12||ap==14||ap==16) nus.push_back(&p);
        }
        if (nus.empty()) return;
        auto* nuRoot = *std::max_element(nus.begin(), nus.end(), [](auto* a, auto* b){ return a->GetEnergy()<b->GetEnergy(); });

        // Start a new tree and add the orginal root, the "hadron blob", the leptons, and photons
        I3MCTreePtr ptree(new I3MCTree);
        using namespace I3MCTreeUtils;
        AddPrimary(*ptree, *nuRoot);
        AppendChild(*ptree, nuRoot->GetID(), *cascade);
        std::set<int> keep = {11,13,15,22};
        for (auto &p : *orig) {
            int ap = std::abs(p.GetPdgEncoding());
            if (ap==12||ap==14||ap==16) {
                if (&p == nuRoot) continue;
                AppendChild(*ptree, nuRoot->GetID(), p);
            }
            else if (keep.count(ap)) {
                I3Particle copy = p;
                if (ap==11||ap==22) {
                    copy.SetShape       (I3Particle::Cascade);
                    copy.SetLocationType(I3Particle::InIce);
                }
                AppendChild(*ptree, nuRoot->GetID(), copy);
            }
        }
        frame->Put("ParamMCTree", ptree);
    }
    void AddHadronTypedMCTree(I3FramePtr frame, const I3MCTreePtr& originalTree) {
        I3MCTreePtr typedTree(new I3MCTree(*originalTree)); // deep copy
        for (auto it = typedTree->begin(); it != typedTree->end(); ++it) {
            const int pdg = std::abs(it->GetPdgEncoding());
            if (pdg == 11 || pdg == 13 || pdg == 15 || // leptons
                pdg == 22 || // photon
                pdg == 12 || pdg == 14 || pdg == 16) // neutrinos
                continue;
            it->SetType(I3Particle::Hadrons);
        }
        frame->Put("HadronsTypedMCTree", typedTree);
    }
    void MakeNCRecoilDarkTree(I3FramePtr frame, const I3MCTreePtr& originalTree) {
        I3MCTreePtr newTree(new I3MCTree(*originalTree));
        
        // Find primary neutrino
        I3Particle* primaryNu = nullptr;
        for (auto it = newTree->begin(); it != newTree->end(); ++it) {
            if (it->GetShape() == I3Particle::Primary && 
                (std::abs(it->GetPdgEncoding()) == 12 || std::abs(it->GetPdgEncoding()) == 14 || std::abs(it->GetPdgEncoding()) == 16)) {
                primaryNu = &(*it);
                break;
            }
        }
        
        if (!primaryNu) {
            frame->Put("NCRecoilDarkTree", newTree);
            return;
        }
        
        // Check children: find charged leptons and recoil neutrinos
        bool hasChargedLepton = false;
        std::vector<I3Particle*> recoilNus;
        
        for (auto it = newTree->begin(); it != newTree->end(); ++it) {
            if (&(*it) == primaryNu || it->GetShape() == I3Particle::Primary) continue;
            
            int pdg = std::abs(it->GetPdgEncoding());
            if (pdg == 11 || pdg == 13 || pdg == 15) hasChargedLepton = true;
            if (pdg == 12 || pdg == 14 || pdg == 16) recoilNus.push_back(&(*it));
        }
        
        // If NC event, darken recoil neutrinos
        if (!hasChargedLepton && !recoilNus.empty()) {
            for (I3Particle* recoil : recoilNus) {
                I3Particle modified = *recoil;
                modified.SetShape(I3Particle::Dark);
                *recoil = modified;
            }
            log_info("NC event: darkened %zu recoil neutrinos", recoilNus.size());
        }
        
        frame->Put("NCRecoilDarkTree", newTree);
    }
}

typedef I3PythiaWrapperModule pythia_wrapper;
namespace { I3_MODULE(pythia_wrapper); }

I3PythiaWrapperModule::~I3PythiaWrapperModule() = default;

// --- Initializing Pythia8 within Icetray, calling the cmnd file ---
void I3PythiaWrapperModule::Configure() {
    GetParameter("PythiaCommands",    pythiaCommands_);
    GetParameter("PythiaCommandPath", pythiaCommandPath_);

    GetParameter("KinematicsKey",     kinematicsKey_);
    GetParameter("VertexKey",         vertexKey_);
    GetParameter("DirectionKey",      directionKey_);
    GetParameter("InteractionKey",    interactionKey_);
    GetParameter("UseFrameKinematics",useFrameKin_);
    GetParameter("YTolerance",        yTolerance_);
    using namespace Pythia8;
    Pythia &pythia = pythia_;
    for (auto const &cmnd : pythiaCommands_) pythia.readString(cmnd);
        if (!pythiaCommandPath_.empty()) {
            LHAPDF6 pdf(&pythia, &pythia.settings, &pythia.logger);
            pythia.readFile(pythiaCommandPath_);
        }
    if (!pythia.init()) {
        log_error("Pythia failed to initialize!");
        return;
    }
    log_info("Pythia initialized!");
    log_info("Wrapper keys: kin=%s vtx=%s dir=%s type=%s",
           kinematicsKey_.c_str(), vertexKey_.c_str(),
           directionKey_.c_str(), interactionKey_.c_str());
}


// --- Generating the DAQ I3MCTree frames by sorting the Pythia8 outputs ---
void I3PythiaWrapperModule::DAQ(I3FramePtr frame) {
    using namespace Pythia8;
    Pythia &pythia = pythia_;

    if (!pythia.next()) {
        log_warn("Pythia failed to generate event.");
        return;
    }
    if (GetIcetrayLogger()->LogLevelForUnit("pythia-wrapper") <= I3LOG_DEBUG)
        pythia.stat();

    // --- prepare containers ---
    I3MCTreePtr mctree(new I3MCTree());
    std::map<int, I3Particle>                      particle_map;
    std::map<int, std::vector<I3MCTree::iterator>> tree_nodes;  // holds *all* clones

    // --- count raw gluons ---> for logging purposes, can be deleted or commented out safely
    int nGluons = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
        if (std::abs(pythia.event[i].id()) == 21) ++nGluons;
    log_info("pythia.event.size() = %d", int(pythia.event.size()));
    log_info("Number of raw Pythia gluons: %d", nGluons);

    // --- stash every Pythia8 particle in an event ---> In Pythia, daughter particles may appear before their parent
    for (int i = 0; i < pythia.event.size(); ++i) {
        auto const &p = pythia.event[i];
        I3Particle tmp;
        tmp.SetEnergy(p.e());
        tmp.SetDir(I3Direction(p.px(), p.py(), p.pz()));
        tmp.SetPos(I3Position(p.xProd()*1e-3, p.yProd()*1e-3, p.zProd()*1e-3));
        tmp.SetTime(p.tProd());
        tmp.SetPdgEncoding(p.id());
        tmp.SetLocationType(I3Particle::InIce);
        //tmp.SetLength(1 * I3Units::m);
        
        //log_info("DBG: part %d → PDG=%d  E=%.6f GeV", i, p.id(), p.e());
        int pdg = std::abs(p.id());
        if ((std::abs(p.id()) == 12 || std::abs(p.id()) == 14 || std::abs(p.id()) == 16) && p.e() >= 9.99e4) {// is neutrino and E >= 100 TeV
            log_info("  → caught neutrino ≥100 TeV");
            tmp.SetShape(I3Particle::Primary);
            tmp.SetLocationType(I3Particle::Anywhere);
        }
        else if (pdg == 13 || pdg == 15) tmp.SetShape(I3Particle::InfiniteTrack);
        else tmp.SetShape(I3Particle::Cascade);
        particle_map[i] = tmp.Clone(); // prevents mixup for Particle ID
    }
    log_info("Stored %zu particles in particle_map", particle_map.size());

    // --- build the tree with duplication ---> assigning the correct order of parents and children
    std::set<int> unresolved;
    for (int i = 0; i < pythia.event.size(); ++i)
        unresolved.insert(i);

    bool progress = true;
    int  attempts = 0, max_attempts = 20;
    while (!unresolved.empty() && progress && attempts++ < max_attempts) {
        progress = false;
        for (auto it = unresolved.begin(); it != unresolved.end();) {
            int i = *it;
            auto const &p = pythia.event[i];
            I3Particle &ip = particle_map[i];

            // skip system (#0)
            if (i == 0) {
                it = unresolved.erase(it);
                continue;
            }

            int m1 = p.mother1(), m2 = p.mother2();
            int st = std::abs(p.status());

            // Case 1: no mother → top‐level
            if (m1 == 0 && m2 == 0) {
                I3Particle child = ip.Clone();            
                auto iter = mctree->insert_after(mctree->end(), child);
                tree_nodes[i].push_back(iter);           
            }
            // Case 2: carbon‐copy
            else if (m1 == m2 && m1 > 0 && tree_nodes.count(m1)) {
                for (auto const &mom : tree_nodes[m1]) {
                    I3Particle child = ip.Clone();         
                    auto iter = mctree->append_child(mom, child);
                    tree_nodes[i].push_back(iter);         
                }
            }
            // Case 3: single mother
            else if (m2 == 0 && m1 > 0 && tree_nodes.count(m1)) {
                for (auto const &mom : tree_nodes[m1]) {
                    I3Particle child = ip.Clone();          
                    auto iter = mctree->append_child(mom, child);
                    tree_nodes[i].push_back(iter);        
                }
            }
            // Case 4 & 5: two mothers
            else if (m2 > m1 && m1 > 0 && tree_nodes.count(m1) && tree_nodes.count(m2)) {                           
                // hang under every m1
                for (auto const &mom1 : tree_nodes[m1]) {
                    I3Particle child = ip.Clone();        
                    auto iter1 = mctree->append_child(mom1, child);
                    tree_nodes[i].push_back(iter1);    
                    if (std::abs(pythia.event[m1].id()) == 12) {
                        log_info("[possible dup‐νe] PDG=%d attached under νₑ (idx=%d), copy 1", p.id(), m1);
                    }
                }
               // and under every m2 (unless the parent is a gluon, skip)
                if (std::abs(p.id()) != 21) {
                    for (auto const &mom2 : tree_nodes[m2]) {
                        I3Particle child2 = ip.Clone();     
                        auto iter2 = mctree->append_child(mom2, child2);
                        tree_nodes[i].push_back(iter2);   
                        if (std::abs(pythia.event[m2].id()) == 12) {
                            log_info("[yes dup‐νe] clone PDG=%d attached under νₑ (idx=%d), copy 2", p.id(), m2);
                        }
                    }
                }
            }
            // Case 6: pick non‐gluon parent
            else if (m1 > m2 && m2 > 0 && tree_nodes.count(m1) && tree_nodes.count(m2)) {
                auto const &p1 = pythia.event[m1];
                int chosen = (std::abs(p1.id()) == 21 ? m2 : m1);
                for (auto const &mom : tree_nodes[chosen]) {
                    I3Particle child = ip.Clone();        
                    auto iter = mctree->append_child(mom, child);
                    tree_nodes[i].push_back(iter);        
                }
            }
            else {
                ++it;
                continue;
            }
        it = unresolved.erase(it);
        progress = true;
        }
    }

    // --- count gluons in tree before pruning ---> For logging purposes, can delete or comment with exception
    size_t gluonsBefore = 0;
    for (auto it = mctree->begin(); it != mctree->end(); ++it) {
        if (std::abs(it->GetPdgEncoding()) == 21)
        ++gluonsBefore;
    }
    log_info("Gluon nodes BEFORE pruning: %zu", gluonsBefore);
    log_info("Total MCTree size before prune: %zu", mctree->size());
    PruneIntermediateGluons(*mctree); // must keep

    size_t gluonsAfter = 0;
    for (auto it = mctree->begin(); it != mctree->end(); ++it) {
        if (std::abs(it->GetPdgEncoding()) == 21)
        ++gluonsAfter;
    }
    log_info("Gluon nodes AFTER pruning:  %zu", gluonsAfter);
    log_info("Total gluon nodes removed: %zu", gluonsBefore - gluonsAfter);
    log_info("Total MCTree size after prune:  %zu", mctree->size());
    RemoveProtonBranch(*mctree); // must keep
    log_info("Total MCTree size after dropping proton: %zu", mctree->size());
    frame->Put("I3MCTree", mctree); ///NEW!!
    MakeHadronizationTree(frame);  ///NEW!!!
    I3MCTreePtr intermediateMCTree(new I3MCTree(*mctree));
    frame->Put("IntermediateMCTree", intermediateMCTree);
    KeepOnlyChildrenParticles(*mctree); // must keep
    log_info("Final MCTree size: %zu", mctree->size());

    // --- offset test ---
    // const I3Position offset(200.0, 100.0, 50.0);
    // for (auto it = mctree->begin(); it != mctree->end(); ++it) {
    //     I3Particle p = *it;
    //      p.SetPos(offset + p.GetPos());
    //      *it = p;
    // }

    // --- push out frame to I3MCTree ---
    frame->Put("I3EventHeader", boost::make_shared<I3EventHeader>());
    frame->Delete("I3MCTree");
    frame->Put("I3MCTree", mctree);
    AddHadronTypedMCTree(frame, mctree); 
    MakeParamCascade(frame);
    MakeParamMCTree(frame);
    MakeNCRecoilDarkTree(frame, mctree);

    PushFrame(frame);

    // --- histogram finals ---> can safely detele or comment out .cxx and .hh
    pdgHistogram_.clear();
    for (auto it = mctree->begin(); it != mctree->end(); ++it) {
        pdgHistogram_[ it->GetPdgEncoding() ]++;
    }
    log_info("Built histogram from %zu tree nodes", mctree->size());
    log_info("[DEBUG] Completed DAQ frame processing");
}

// --- creating a csv for particle count ---> can safely delete or comment out in .cxx and .hh
void I3PythiaWrapperModule::Finish() {
    std::ofstream out("particle_counts.csv");
    out << "PDG_ID,Count\n";
    for (auto const & [pdg, count] : pdgHistogram_)
        out << pdg << "," << count << "\n";
    out.close();
    log_info("Saved particle histogram to particle_counts.csv");
}

bool I3PythiaWrapperModule::ShouldDoProcess(I3FramePtr) { return true; }

