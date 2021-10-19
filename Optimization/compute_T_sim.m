function T_sim = compute_T_sim(FT,aT,MAinterp,Misc)

                T_sim = (FT'*reshape(MAinterp,Misc.NMuscles,Misc.nDOF))' + Misc.Topt*aT;
            

end