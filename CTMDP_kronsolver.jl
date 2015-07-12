using kronfun
using HDF5, JLD

using CTMDP


#It is assumed that someone else is providing the following:
#g_allStates: an array of symbols representing all possible substates that we can be in
#g_sn: a dictionary of type (Symbol => Int64) which goes from substate to index
#legalActions(S): function which tells us which actions are allowed!
#g_nullAct: the noaction case
