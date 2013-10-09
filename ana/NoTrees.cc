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


//} // namespace contrib
//
//FASTJET_END_NAMESPACE
