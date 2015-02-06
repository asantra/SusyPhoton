{
 Double_t ratioWeight = (rewByAFloat) ? ratio->GetBinContent(ratio->FindBin(rewVar)) : ratio->GetBinContent(ratio->FindBin(rewVar_int));
Double_t ratioWeightErr = (rewByAFloat) ? ratio->GetBinError(ratio->FindBin(rewVar)) : ratio->GetBinError(ratio->FindBin(rewVar_int));

if(ratioWeight == 0.) continue;

Double_t oldError = (outIsAFloat) ? rew->GetBinError(rew->FindBin(outVar)) : rew->GetBinError(rew->FindBin(outVar_int));

if(outIsAFloat) rew->Fill(outVar, 1.0*ratioWeight);
else rew->Fill(outVar_int, 1.0*ratioWeight);

Double_t newError = ratioWeightErr*ratioWeightErr / (ratioWeight*ratioWeight);

if(outIsAFloat) rew->SetBinError(rew->FindBin(outVar), sqrt(oldError*oldError + newError*newError));
else rew->SetBinError(rew->FindBin(outVar_int), sqrt(oldError*oldError + newError*newError));
}
