"""
A set of functions for the estimation of the methylation and acetylation factors 
                as outlined in the Pratx The EMBO Journal e113595 You can use this
                for a single gene or multiple genes if the genes expression or the 
                immunoprecipitation is supplied in a csv file. It uses a vector 
                assignment by using the dataframe approach as compared to the 
                direct indiviual math implementation for faster assignment and 
                processing.
plotting functions are also incorporated. The variables should have a single column 
in the csv file otherwise change the code using the column assignment [[]]. 
"""

def EstMeRetained(eGFP,eGFPMeth):
    """
    python implementation of the histone turnover estimation. 
    You can read a file and analyze as many genes as you want.
    """
    import pandas as pd
    GFP_meth24expression_read = pd.read_csv(eGFPMeth,sep = "\t")
    GFP_expression_read = pd.read_csv(eGFP, sep = "\t")
    MeRetained = pd.DataFrame(GFP_meth24expression_read/GFP_expression_read)
    return MeRetained.to_list()

def alpha(GFP_t_start,GFP_t_end, 
                       H3_t_start, 
                              H3_t_end, plot="false"):
    """
    GFP_t_start and GFP_t_end are the start and the end
    timepoints of the GFP tagged.
    H3_t_start and H3_t_end are the start and the end
    timepoints of the H3 untagged.
    """
    if plot == "false":
        import pandas as pd
        GFP_t_start = pd.read_csv(GFP_t_start, sep = "\t")
        GFP_t_end = pd.read_csv(GFP_t_end, sep = "\t")
        H3_t_start = pd.read_csv(H3_t_start, sep = "\t")
        H3_t_end = pd.read_csv(H3_t_end, sep = "\t")
        end_point = GFP_t_end/H3_t_end 
        start_point = GFP_t_start/H3_t_start
        calculated_alpha = (1-GFP_t_start/H3_t_start)*end_point/start_point
        return calculated_alpha
    elif plot == "bar":
        import pandas as pd
        GFP_t_start = pd.read_csv(GFP_t_start, sep = "\t")
        GFP_t_end = pd.read_csv(GFP_t_end, sep = "\t")
        H3_t_start = pd.read_csv(H3_t_start, sep = "\t")
        H3_t_end = pd.read_csv(H3_t_end, sep = "\t")
        end_point = GFP_t_end/H3_t_end 
        start_point = GFP_t_start/H3_t_start
        calculated_alpha = (1-GFP_t_start/H3_t_start)*end_point/start_point
        return pd.DataFrame(calculated_alpha, columns = ["alpha"]).plot.bar()
    else:
        print("no option selected")

def Met0(H3_start, H3_meth):
    """
    """
    import pandas as pd
    start = pd.read_csv(H3_start, sep = "\t")
    end = pd.read_csv(H3_meth, sep = "\t")
    return pd.DataFrame(start/end, columns = ["Met0"]).to_list() 
