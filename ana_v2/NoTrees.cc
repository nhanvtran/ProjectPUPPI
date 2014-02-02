// NoTrees Package
// Questions/Comments? danbert@mit.edu jthaler@mit.edu
//
// Copyright (c) 2013
// Daniele Bertolini, and Jesse Thaler
//
// $Id: $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "NoTrees.hh"
#include <sstream>

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//namespace contrib{


//----------------------------------------------------------------------
//Helper for SelectorEventTrimmer. 
//Class derived from SelectorWorker: pass, terminator, applies_jet_by_jet, and description are overloaded. 
//It can't be applied to an individual jet since it requires information about the neighbourhood of the jet.


class SW_EventTrimmer: public SelectorWorker {
public:
   
   SW_EventTrimmer(double Rjet, double ptcut, double Rsub, double fcut) : _Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut) {};

   virtual bool pass(const PseudoJet &) const {
	if (!applies_jet_by_jet())	
     	throw Error("Cannot apply this selector worker to an individual jet");
	return false;
   }

   virtual void terminator(vector<const PseudoJet *> & jets) const {
	//copy pointers content to a vector of PseudoJets. Check if the pointer has been already nullified. 
	vector<PseudoJet> my_jets;
	vector<unsigned int> indices; 
	for(unsigned int i=0; i < jets.size(); i++){
		if(jets[i]){
			indices.push_back(i);
			my_jets.push_back(*jets[i]);
  		} 	
	}
	
   	for(unsigned int i=0; i < my_jets.size(); i++){
 		double pt_Rjet = pt_within_R(my_jets,my_jets[i],_Rjet);
		double pt_Rsub = pt_within_R(my_jets,my_jets[i],_Rsub);   
		if(pt_Rjet <= _ptcut || pt_Rsub/pt_Rjet <= _fcut) jets[indices[i]] = NULL;
	}
   }

   virtual bool applies_jet_by_jet() const {return false;}

   virtual string description() const {
	ostringstream ostr;
	ostr << "Event trimmer with R_jet= " << _Rjet <<", pT_cut=" << _ptcut;
        ostr << ", R_sub=" << _Rsub << ", and f_cut=" << _fcut;
	return ostr.str();
   }

private:

  double _Rjet, _ptcut, _Rsub, _fcut;

}; 


Selector SelectorEventTrimmer(double Rjet, double ptcut, double Rsub, double fcut){
   return Selector(new SW_EventTrimmer(Rjet,ptcut,Rsub,fcut));
} 



//----------------------------------------------------------------------

double JetMultiplicity::result(const vector<PseudoJet> & particles) const {
	
   double answer = 0.0;
  
   if(_trim == true) {
   
   	for (unsigned int i = 0; i < particles.size(); i++) {
		
		double pt_Rjet = pt_within_R(particles,particles[i],_Rjet);
		double pt_Rsub = pt_within_R(particles,particles[i],_Rsub);
	
		if(pt_Rjet >= _ptcut && pt_Rsub/pt_Rjet >= _fcut) answer += particles[i].pt()/pt_Rjet;
	}
   }
   
   else if(_trim == false) {
   
   	for (unsigned int i = 0; i < particles.size(); i++) {
		
		double pt_Rjet = pt_within_R(particles,particles[i],_Rjet);
		if(pt_Rjet >= _ptcut) answer += particles[i].pt()/pt_Rjet;
	}
   }

   return(answer);
}


double TransverseEnergy::result(const vector<PseudoJet> & particles) const {

   double answer = 0.0;
  
   if(_trim == true) {
   
   	for (unsigned int i = 0; i < particles.size(); i++) {
		
		double pt_Rjet = pt_within_R(particles,particles[i],_Rjet);
		double pt_Rsub = pt_within_R(particles,particles[i],_Rsub);
	
		if(pt_Rjet >= _ptcut && pt_Rsub/pt_Rjet >= _fcut) answer += particles[i].pt();
	}
   }
   
   else if(_trim == false) {
   
   	for (unsigned int i = 0; i < particles.size(); i++) {
		
		double pt_Rjet = pt_within_R(particles,particles[i],_Rjet);
		if(pt_Rjet >= _ptcut) answer += particles[i].pt();
	}
   }

   return(answer);
}


double MissingTransverseEnergy::result(const vector<PseudoJet> & particles) const {

   PseudoJet p_tot(0,0,0,0);
   
   if(_trim == true) {
   
   	for (unsigned int i = 0; i < particles.size(); i++) {
		
		double pt_Rjet = pt_within_R(particles,particles[i],_Rjet);
		double pt_Rsub = pt_within_R(particles,particles[i],_Rsub);
	
		if(pt_Rjet >= _ptcut && pt_Rsub/pt_Rjet >= _fcut) p_tot += particles[i];
	}
   }
   
   else if(_trim == false) {
   
   	for (unsigned int i = 0; i < particles.size(); i++) {
		
		double pt_Rjet = pt_within_R(particles,particles[i],_Rjet);
		if(pt_Rjet >= _ptcut) p_tot += particles[i];
	}
   }

   return(p_tot.pt());
}



string JetMultiplicity::description() const {
  ostringstream oss;
  if(_trim == false) oss << "Treeless Jet multiplicity with R_jet=" << _Rjet << ", and pT_cut=" << _ptcut;
  else if(_trim == true) {
  	oss << "Treeless Jet multiplicity with R_jet=" << _Rjet << ", and pT_cut=" << _ptcut;
	oss << " + Trimming with R_sub=" << _Rsub <<", and fcut=" << _fcut;
  }
  return oss.str();
}

string TransverseEnergy::description() const {
  ostringstream oss;
  if(_trim == false) oss << "Treeless Transverse Energy with R_jet=" << _Rjet << ", and pT_cut=" << _ptcut;
  else if(_trim == true) {
  	oss << "Treeless Transverse Energy with R_jet=" << _Rjet << ", and pT_cut=" << _ptcut;
	oss << " + Trimming with R_sub=" << _Rsub <<", and fcut=" << _fcut;
  }
  return oss.str();
}

string MissingTransverseEnergy::description() const {
  ostringstream oss;
  if(_trim == false) oss << "Treeless Missing Transverse Energy with R_jet=" << _Rjet << ", and pT_cut=" << _ptcut;
  else if(_trim == true) {
  	oss << "Treeless Missing Transverse Energy with R_jet=" << _Rjet << ", and pT_cut=" << _ptcut;
	oss << " + Trimming with R_sub=" << _Rsub <<", and fcut=" << _fcut;
  }
  return oss.str();
}


//----------------------------------------------------------------------

double pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R){

   fastjet::Selector sel = fastjet::SelectorCircle(R);
   sel.set_reference(centre);
   vector<PseudoJet> near_particles = sel(particles);
   double answer = 0.0;

   for(unsigned int i=0; i<near_particles.size(); i++){
	answer += near_particles[i].pt();
   }
   return(answer);
}

PseudoJet flow_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R){

   fastjet::Selector sel = fastjet::SelectorCircle(R);
   sel.set_reference(centre);
   vector<PseudoJet> near_particles = sel(particles);
   PseudoJet flow;
   for(unsigned int i=0; i<near_particles.size(); i++){
     //if(flow.pt() < near_particles[i].pt()) flow = near_particles[i];
     flow += near_particles[i];
   }
   return flow;
}

double var_within_R(int iId, const vector<PseudoJet> & particles, const PseudoJet& centre, double R){
  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  vector<PseudoJet> near_particles = sel(particles);
  double var = 0;
  double lSumPt = 0;
  if(iId == 5 || iId == 6) for(unsigned int i=0; i<near_particles.size(); i++) lSumPt += near_particles[i].pt();
  for(unsigned int i=0; i<near_particles.size(); i++){
    double pDEta = near_particles[i].eta()-centre.eta();
    double pDPhi = fabs(near_particles[i].phi()-centre.phi());
    if(pDPhi > 2.*3.14159265-pDPhi) pDPhi =  2.*3.14159265-pDPhi;
    double pDR = sqrt(pDEta*pDEta+pDPhi*pDPhi);
    if(pDR  < 0.001) continue;
    if(pDR  <  0.05) pDR = 0.05;
    if(pDR == 0) continue;
    if(iId == 0) var += pDR;
    if(iId == 1) var += near_particles[i].pt();
    if(iId == 2) var += log(pDR*near_particles[i].pt());
    if(iId == 3) var += log(near_particles[i].pt()/sqrt(pDR));
    if(iId == 4) var += log(near_particles[i].pt()/pDR);
    if(iId == 5) var += log(near_particles[i].pt()/pDR/lSumPt);
    if(iId == 6) var += log(near_particles[i].pt()/pDR);
    if(iId == 7) var += log(near_particles[i].pt()/pDR);
  }
  return var;
}

double puppi_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R, bool iMass){

    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    
//    // pT_{i,R_puppi}
//    double ptIR = 0.0;
//    for(unsigned int i=0; i<near_particles.size(); i++){
//        ptIR += near_particles[i].pt();
//    }
    
    double puppi = 0;
    for(unsigned int i=0; i<near_particles.size(); i++){
        double pDEta = near_particles[i].eta()-centre.eta();
        double pDPhi = fabs(near_particles[i].phi()-centre.phi());
        if(pDPhi > 2.*3.14159265-pDPhi) pDPhi =  2.*3.14159265-pDPhi;
        double pDR = sqrt(pDEta*pDEta+pDPhi*pDPhi);
        if(pDR  < 0.001) continue;
        if(pDR  <  0.05) pDR = 0.05;
        if(pDR == 0) continue;
        double lE =  near_particles[i].E() + centre.E();
        //if(lE < 0.1) continue;
        double lM =  (near_particles[i] + centre).m();
        double lP =  sqrt((near_particles[i] + centre).pt2()+((near_particles[i] + centre).pz()*(near_particles[i] + centre).pz()));
        double lDot = near_particles[i].E()*centre.E() - near_particles[i].px()*centre.px()-near_particles[i].py()*centre.py()-near_particles[i].pz()*centre.pz();
        PseudoJet lSum = centre + near_particles[i];
        double lFSum = near_particles[i].E()*lSum.E() - near_particles[i].px()*lSum.px()-near_particles[i].py()*lSum.py()-near_particles[i].pz()*lSum.pz();
        // cout << "---> " << lDot << " -- " << lFSum << endl;
        lDot = max(lFSum,lDot);
        if(lE < lP) lM = 0.01; //cout << "===> Problem " <<  pDR << " -- " << lE << " -- " << lP << endl;
        if(lE > lP) lM = sqrt(lE*lE-lP*lP);
        if(lM < 0.01) lM = 0.01;
        
        //double puppi_wExp1 = 2.*log(max(near_particles[i].pt(),centre.pt())/pDR);
        //double puppi_wExp2 = 2.*log(max(near_particles[i].pt(),centre.pt())*max(near_particles[i].pt(),centre.pt())/pDR/pDR);        
        //std::cout << "puppi_wExp1 = " << puppi_wExp1 << ", puppi_wExp2 = " << puppi_wExp2 << std::endl;
        
        //=>if(!iMass) puppi += min(near_particles[i].pt(),centre.pt())*pDR;//log(pDR/lE);
        //if(!iMass) puppi += log(max(near_particles[i].pt(),centre.pt())/pDR);//log(pDR/lE);
       
        //puppi += log(near_particles[i].pt()/pDR);
        puppi += log(near_particles[i].pt()*centre.pt()/pDR/pDR);
        
        if(iMass) puppi += max(near_particles[i].pt(),centre.pt())*max(near_particles[i].pt(),centre.pt())/pDR/pDR;//*(near_particles[i]+centre).pt()*(near_particles[i]+centre).pt();
        
    }
    
//    double weight = centre.pt()/ptIR;
//    std::cout << "weight = " << weight << ", near_particles.size() = " << near_particles.size() << std::endl;
//    puppi *= centre.pt()/ptIR;
    
    return puppi;
}

//} // namespace contrib
//
//FASTJET_END_NAMESPACE
