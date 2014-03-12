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

#ifndef __FASTJET_CONTRIB_NOTREES_HH__
#define __FASTJET_CONTRIB_NOTREES_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"

#include <string>
#include <algorithm>

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
using namespace fastjet;


//namespace contrib{

//----------------------------------------------------------------------
//Event trimmer. It is implemented as a selector.

Selector SelectorEventTrimmer(double Rjet, double ptcut, double Rsub, double fcut);

//----------------------------------------------------------------------
//In current version of FastJet there is no standard interface for a function that takes a vector of PseudoJets as argument. 
//This class serves as a base class and it is the analog of FunctionOfPseudoJet, it only differs in taking a vector of PseudoJets as argument.
//Should a standard interface become available, this class can be removed and the classes below
//(JetMultiplicity,TransverseEnergy,MissingTransverseEnergy) should be derived from the standard base class. 

template<typename TOut>
class MyFunctionOfPseudoJets{
public:

   MyFunctionOfPseudoJets(){}

   MyFunctionOfPseudoJets(const TOut &constant_value);
 
   virtual ~MyFunctionOfPseudoJets(){}

   virtual std::string description() const{ return "";}

   virtual TOut result(const std::vector<PseudoJet> &pjs) const = 0;

   TOut operator()(const std::vector<PseudoJet> &pjs) const { return result(pjs);}
};


//----------------------------------------------------------------------
// JetMultiplicity defines treeless jet/subjet counters. Optional built-in trimmer.
class JetMultiplicity: public MyFunctionOfPseudoJets<double> {

public:

   JetMultiplicity(){};
  
   JetMultiplicity(double Rjet, double ptcut) : _Rjet(Rjet), _ptcut(ptcut), _trim(false) {};
   
   JetMultiplicity(double Rjet, double ptcut, double Rsub, double fcut) : 
	_Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _trim(true) {};

   virtual ~JetMultiplicity(){}

   double result(const std::vector<PseudoJet> & particles) const;
  
   std::string description() const;

    
private:
 
   double _Rjet, _ptcut, _Rsub, _fcut;
   bool _trim;

};


//----------------------------------------------------------------------
// TransverseEnergy defines treeless scalar pT sum. Optional built-in trimmer.

class TransverseEnergy: public MyFunctionOfPseudoJets<double> {

public:
  
   TransverseEnergy(double Rjet, double ptcut) :  _Rjet(Rjet), _ptcut(ptcut), _trim(false) {};
   TransverseEnergy(double Rjet, double ptcut, double Rsub, double fcut) : 
	_Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _trim(true) {};

   virtual ~TransverseEnergy(){}

   double result(const std::vector<PseudoJet> & particles) const;
 
   std::string description() const;

   
private:
 
   double _Rjet, _ptcut, _Rsub, _fcut;
   bool _trim;
   
};


//----------------------------------------------------------------------
// MissingTransverseEnergy defines treeless pT of the vector sum of jets momenta. Optional built-in trimmer.

class MissingTransverseEnergy: public MyFunctionOfPseudoJets<double> {

public:
  
   MissingTransverseEnergy(double Rjet, double ptcut) :_Rjet(Rjet), _ptcut(ptcut), _trim(false) {};
   
   MissingTransverseEnergy(double Rjet, double ptcut, double Rsub, double fcut) : 
   	_Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _trim(true) {};

   virtual ~MissingTransverseEnergy(){}

   double result(const std::vector<PseudoJet> & particles) const;
 
   std::string description() const;

   
private:
 
   double _Rjet, _ptcut, _Rsub, _fcut;
   bool _trim;

};


//----------------------------------------------------------------------
//Function pt_within_R sums pT of a subset of particles. The subset is defined by a circle of radius R around centre.

double var_within_R(int iId, const vector<PseudoJet> & particles, const PseudoJet& centre, double R);
double pt_within_R(const std::vector<PseudoJet> & particles, const PseudoJet& centre, double R);
PseudoJet flow_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R);
double puppi_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R,bool iMass);


//} // namespace contrib

//FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NOTREES_HH__
