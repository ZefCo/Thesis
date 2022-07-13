from importlib.resources import path
import pandas
import pathlib
import numpy
import re
import random
# import openpyxl
import itertools
import sequence_alignments as sa
from pca import pca
import plotly.express as px
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
# from rpy2.robjects import 
from rpy2.robjects.conversion import localconverter



class EE:
    def __init__(self) -> None:

        self.length = 10
        self.master_freq = ["33P", "33Ran", "55P", "55Ran"]

        self.FFP_stat, self.FFR_stat, self.TTP_stat, self.TTR_stat = pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame()
        self.Fre_xlsx, self.FFP_xlsx, self.FFR_xlsx, self.TTP_xlsx, self.TTR_xlsx, self.Fre_xlsx, self.Stat_xlsx = {}, {}, {}, {}, {}, {}, {}
        
        self.weights()
        self.import_settings()
        self.import_data()
        self.exon2exon()
        self.pca_output()
        # PCA Analysis of Mean Median Mode etc
            # How to handel Mode? Should I? For now don't bother.
            # Mean, Median, Max, Min, Var, Skew, Kutosis
            # Mode will have to be handeled carefully
        # Print Covar matrix
        # Print PC**3
        self.lda_output()
        self.write2excel()


    def convert_2_R(self, dataframe: pandas.DataFrame):
        '''
        '''
        # So this works but I have to think of a more general approach at some point: everything is a numpy.float64 to start and I need to convert that to
        # a float to get it into R.
        try:
            with localconverter(robjects.default_converter + pandas2ri.converter):
                dataframe = robjects.conversion.py2rpy(dataframe)
        except AttributeError:
            columns = dataframe.columns
            for col in columns:
                # print(dataframe[col].dtype)
                if dataframe[col].dtype == 'category':
                    # Idk... I wrote this as a print statement cause I thought maybe I had to do more... but I don't so... whatever
                    pass
                else:
                    dataframe[col] = dataframe[col].astype(float)
            
            with localconverter(robjects.default_converter + pandas2ri.converter): 
                dataframe = robjects.conversion.py2rpy(dataframe)
            
        except Exception as e:
            print(type(e))


        return dataframe




    def frequency_results(self, input_dataframe: pandas.DataFrame, name: str):
        unique_gnum = set()
        colnames = input_dataframe.columns

        for col in colnames:
            unique_columns = set(input_dataframe[col])
            unique_gnum = unique_gnum.union(unique_columns)
        unique_gnum = sorted(unique_gnum)

        unique_gnum = pandas.Series(0, index = unique_gnum)
        # unique_stat = pandas.Series(0, index = ["Mean", "Median", "Mode", "Range", "SD", "Var"])

        for col in colnames:
            unique_columns = input_dataframe[col].value_counts()
            unique_gnum = unique_gnum.add(unique_columns, fill_value = 0)

        unique_gnum.name = name
        
        # print(unique_gnum)

        return unique_gnum




    def exon2exon(self):
        '''
        '''
        # print(self.utdata)
        # print(self.exon_refseq)
        fg_row, _ = self.utdata.shape

        for frow in range(fg_row):
            fusion_of_interest = self.utdata.iloc[frow, :]
            # What comes in is a string, but it looks like a Python list: so it has to be stripped of the thing appearence of the list, then split on the ,
            # to be turned into an actual list
            head_gene, tail_gene = fusion_of_interest['Hgene'], fusion_of_interest['Tgene']
            unique_identifier = f"{head_gene}_{tail_gene}"
            
            # hnames: str = re.sub(r'\[|\]|\'', '', fusion_of_interest['HNames'])
            hnames: list = re.split(',', re.sub(self.junk, '', fusion_of_interest['HNames']))
            # tnames: str = re.sub(r'\[|\]|\'', '', fusion_of_interest['TNames'])
            tnames: list = re.split(',', re.sub(self.junk, '', fusion_of_interest['TNames']))
            # Add in a check to see if they are the same. If they are, skip, if not, keep going

            if len(set(hnames).intersection(set(tnames))) > 0:
                print(f"Check {unique_identifier} for issues: inserction is not zero")
                continue
                
            head_exon_subseq = self.exon_refseq[self.exon_refseq['name'].isin(hnames)]
            tail_exon_subseq = self.exon_refseq[self.exon_refseq['name'].isin(tnames)]

            FFP_frame, FFR_rando, FFP_prime, TTP_frame, TTR_rando, TTP_prime = pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame()

            for hrow in range(head_exon_subseq.shape[0]):
                his_name = head_exon_subseq.iloc[hrow, :]['name']
    #             # Get the head exon count for each isoform
                hexon_count = head_exon_subseq.iloc[hrow, :]['exonCount']
                
    #             # Generate a name list from this row
                head_gene_names = [f'{head_gene}_{his_name}_exon{e}' for e in range(1, hexon_count + 1)]
    #             # print(f"\t{len(head_gene_names)}")

                # This seems backwards but draw it out: the 5' is the last few nucelotides of the sequence, while the 3' is the first few.
                # Trust me: it makes sense when you draw it out.
                # The reason that it's taking a row of the dataframe then taking the column names is because I don't know what numerical index
                # each column has. I grab the whole row, then figure out which columns I need based on the name. You could just count from the index column over
                # to figure out which column is needed, but this was faster to write.
                hexon_5prime = [head_exon_subseq.iloc[hrow, :][f'Exon_{e}'][0 : self.length] for e in range(1, hexon_count + 1)]
                hexon_3prime = [head_exon_subseq.iloc[hrow, :][f'Exon_{e}'][len(head_exon_subseq.iloc[hrow, :][f'Exon_{e}']) - self.length : len(head_exon_subseq.iloc[hrow, :][f'Exon_{e}'])][::-1] for e in range(1, hexon_count + 1)]
                # Sorry about the last line, probably shouldn't have used list comprehension. It's grabbing the end of the string, and the length of the string varries from exon to exon, and then reverses it. It's really just:
                # [grab row at [location][length of string - overall length: to : end of string][reversed] for every exon in row]

                for trow in range(tail_exon_subseq.shape[0]):
                    tis_name = tail_exon_subseq.iloc[trow, :]['name']
                    TTR_subrando, TTP_subframe, FFP_subframe, FFR_subrando = pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame()

                #     # Get the tail exon count
                    texon_count = tail_exon_subseq.iloc[trow, :]['exonCount']
                    
                #     # Generate a name list
                    tail_gene_names = [f'{tail_gene}_{tis_name}_exon{e}' for e in range(1, texon_count + 1)]
                #     # old_tail_gene_names = tail_gene_names
                #     # print(f"\t\t{len(tail_gene_names)}")

                    # Getting all sequences in a list for 5' and 3'
                    texon_5prime = [tail_exon_subseq.iloc[trow, :][f'Exon_{e}'][0 : self.length] for e in range(1, texon_count + 1)]
                    texon_3prime = [tail_exon_subseq.iloc[trow, :][f'Exon_{e}'][len(tail_exon_subseq.iloc[trow, :][f'Exon_{e}']) - self.length : len(tail_exon_subseq.iloc[trow, :][f'Exon_{e}'])][::-1] for e in range(1, texon_count + 1)]
                    # print(f"\t\t{len(texon_5prime)}, {len(texon_3prime)}")

                    for i, hexon5p in enumerate(hexon_5prime):
                        skipping = 0
                        new_d5prime, new_d5rando, new_d3prime, new_d3rando = [], [], [], []

                        hexon3p = hexon_3prime[i]

                        for j, texon5p in enumerate(texon_5prime):

                            texon3p = texon_3prime[j]

                            # Original way of dealing with different lengths: do not align them and just set the values
                            # to 1. Makes sense, it's quick, but it's also dirty. Avoids having length and index errors in sequences.
                            if (len(hexon5p) != self.length) or (len(texon5p) != self.length):
                                skipping += 1
                                # print(f'hexon length: {len(hexon5p)}\ttexon length: {len(texon5p)}')
                                # # print(f"Something weird at index {i} jndex {j}")
                                # print("One of the above sequences is not of the appropriate length\nSetting G number to 1 and moving on")
                                # print("Consider using a modified alignment sequnce for scoring")

                                g_5prime, g_5rando, _ = 1, 1, 1
                                g_3prime, g_3rando, _ = 1, 1, 1

                            else:
                                g_5prime, g_5rando, _ = self.random_comparison(texon5p, hexon5p)
                                g_3prime, g_3rando, _ = self.random_comparison(hexon3p, texon3p)

                            new_d5prime.append(g_5prime), new_d5rando.append(g_5rando), new_d3prime.append(g_3prime), new_d3rando.append(g_3rando)


                        new_5prime_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d5prime)), index = pandas.Index([head_gene_names[i]]))
                        FFP_subframe = pandas.concat([FFP_subframe, new_5prime_row], axis = 0)

                        new_5rando_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d5rando)), index = pandas.Index([head_gene_names[i]]))
                        FFR_subrando = pandas.concat([FFR_subrando, new_5rando_row], axis = 0)

                        new_3prime_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d3prime)), index = pandas.Index([head_gene_names[i]]))
                        TTP_subframe = pandas.concat([TTP_subframe, new_3prime_row], axis = 0)

                        new_3rando_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d3rando)), index = pandas.Index([head_gene_names[i]]))
                        TTR_subrando = pandas.concat([TTR_subrando, new_3rando_row], axis = 0)
                        
                        if skipping > 0:
                            print(f"Skipped a total of {skipping} Exon comparisions - {unique_identifier}: {his_name}_{tis_name} - length was not appropriate")

                    FFP_frame = pandas.concat([FFP_frame.stack(), FFP_subframe.stack()], axis = 0).unstack()
                    TTP_frame = pandas.concat([TTP_frame.stack(), TTP_subframe.stack()], axis = 0).unstack()
                    FFR_rando = pandas.concat([FFR_rando.stack(), FFR_subrando.stack()], axis = 0).unstack()
                    TTR_rando = pandas.concat([TTR_rando.stack(), TTR_subrando.stack()], axis = 0).unstack()
            
            # Taking the transpose of the FFP porition, because the Tail is the second gene and I want that one as the row labels
            FFP_frame = FFP_frame.T
            FFR_rando = FFR_rando.T

            self.FFP_xlsx[unique_identifier] = FFP_frame
            self.TTP_xlsx[unique_identifier] = TTP_frame
            self.FFR_xlsx[f"{unique_identifier}_Random"] = FFR_rando
            self.TTR_xlsx[f"{unique_identifier}_Random"] = TTR_rando

            ffp_substat = self.meta_stats(FFP_frame, f"{unique_identifier}_55P")
            ffr_substat = self.meta_stats(FFR_rando, f"{unique_identifier}_55Ran")
            ttp_substat = self.meta_stats(TTP_frame, f"{unique_identifier}_33P")
            ttr_substat = self.meta_stats(TTR_rando, f"{unique_identifier}_33Ran")

            self.FFP_stat = pandas.concat([self.FFP_stat, ffp_substat.to_frame().T], axis = 0)
            self.FFR_stat = pandas.concat([self.FFR_stat, ffr_substat.to_frame().T], axis = 0)
            self.TTP_stat = pandas.concat([self.TTP_stat, ttp_substat.to_frame().T], axis = 0)
            self.TTR_stat = pandas.concat([self.TTR_stat, ttr_substat.to_frame().T], axis = 0)

            # stat_frame = pandas.concat([FFP_stat.to_frame().stack(),
            #                             FFR_stat.to_frame().stack(),
            #                             TTP_stat.to_frame().stack(),
            #                             TTR_stat.to_frame().stack()], axis = 1).unstack()


            FFP_freq = self.frequency_results(FFP_frame, f"{unique_identifier}_55P")
            FFR_freq = self.frequency_results(FFR_rando, f"{unique_identifier}_55Ran")
            TTP_freq = self.frequency_results(TTP_frame, f"{unique_identifier}_33P")
            TTR_freq = self.frequency_results(TTR_rando, f"{unique_identifier}_33Ran")

            # So this is important: the concatation order is to output the final columns ALPHABETICALLY, which means it's going to be 
            # [UniqueID]_33P    [UniqueID]_33ran    [UniqueID]_55P  [UniqueID]_33Ran
            # So even though I'm putting them in a different order they will output in this order. That's the order that has to be preserved when doing the
            # running totals.
            frequency_frame = pandas.concat([FFP_freq.to_frame().stack(), 
                                                FFR_freq.to_frame().stack(), 
                                                TTP_freq.to_frame().stack(), 
                                                TTR_freq.to_frame().stack()], axis = 0).unstack()
            # print(frequency_frame)
            # frequency_meta = pandas.concat([FFP_stat.to_frame().stack(),
            #                                 FFR_stat.to_frame().stack(),
            #                                 TTP_stat.to_frame().stack(),
            #                                 TTR_stat.to_frame().stack()], axis = 0).unstack()
            # print(frequency_meta)

            self.Fre_xlsx[unique_identifier] = frequency_frame

            frequency_add = frequency_frame.fillna(0)
            frequency_add.columns = self.master_freq
            self.Frequency = self.Frequency.add(frequency_add, fill_value = 0)
            # print(Frequency)

            self.Stat_xlsx['T5P_H5P'], self.Stat_xlsx['T5R_H5R'], self.Stat_xlsx['H3P_T3P'], self.Stat_xlsx['H3R_T3R'] = self.FFP_stat, self.FFR_stat, self.TTP_stat, self.TTR_stat

            # There's probably a much smarter way to handel this. This chunk of code was cut and pasted from the self.pca_output() because I realized I needed to handle it twice, once for
            # pca and once for lda, so I just moved it here. But not I have this setnames which is not handeled properly... Does R convert a catagory to a factor? I may have to adjust the
            # R2P function later
            self.mmm_frame = pandas.DataFrame()
            self.setnames = []
            for sheetname, sheet in self.Stat_xlsx.items():
                self.mmm_frame = pandas.concat([self.mmm_frame, sheet], axis = 0)

                tempnames = [sheetname for _ in range(sheet.shape[0])]
                self.setnames = self.setnames + tempnames

            self.mmm_frame = self.mmm_frame[["Mean", "Median", "Max", "Min", "Var", "Skew", "Kurtosis"]]



    def g_number(self, agene, bgene):
        '''
        '''
        gdelta = numpy.array([0 if na == bgene[i] else 1 for i, na in enumerate(agene)])

        g = numpy.sum(self.weight * gdelta)

        return g



    def import_data(self):
        '''
        '''
        with open(pathlib.Path.cwd() / "Data_Files" / f"{self.ut_file_name}") as utdata_file:
            self.utdata = pandas.read_csv(utdata_file, sep = ',', header = 0)

        with open(pathlib.Path.cwd() / "Data_Files" / f"{self.genome_file}") as genomeseq_file:
            self.exon_refseq = pandas.read_csv(genomeseq_file, sep = ',', header=0)

        

    def import_settings(self):
        '''
        '''
        self.junk = r'\[|\]|\'|\s'
        self.outputpath = pathlib.Path.cwd() / "Data_Files" / "GNum" / f"EEv4_Len_{self.length}"
        self.genome = 'hg19'

        try:
            self.outputpath.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            print("Folder is already here")
        else:
            print(f"Path {self.outputpath} was created")

        self.ut_file_name = 'New_UTData_cds.csv'
        self.genome_file = 'Blast_Genes_HG19.csv'

        self.write_files: bool = True
        self.pca_output_on: bool = True
        self.lca_output_on: bool = True

        # print(f"Head Info\t{head_gene}\t{head_chr}\t{head_strand}")
        # print(F"Tail Info\t{tail_gene}\t{tail_chr}\t{tail_strand}")
        # print("####")

    
    def lda_output(self):
        '''
        '''
        lda_script = '''
        library(MASS)
        function(lataframe) {
            model <- lda(Set ~ ., data = lataframe)
            plotting_data <- predict(model)
            plotting_data <- data.frame(x = plotting_data$x[, 1], y = plotting_data$x[, 2], z = plotting_data$x[, 3], Set = plotting_data$class)

            return(plotting_data)
        }
        '''
        lda = robjects.r(lda_script)

        self.mmm_frame["Set"] = pandas.Categorical(self.setnames)
        mmm_frame = self.convert_2_R(self.mmm_frame)

        mmm_lda = lda(mmm_frame)

        with localconverter(robjects.default_converter + pandas2ri.converter):
            mmm_lda: pandas.DataFrame = robjects.conversion.rpy2py(mmm_lda)
        mmm_lda["Fusion"] = self.mmm_frame.index


        lda_plot = px.scatter(mmm_lda, x = 'x', y = 'y', color = 'Set', symbol = 'Set', title = f"LDA of MMM", hover_data = {"Fusion": True})
        lda_plot.show()
        lda_plot.write_html(pathlib.Path.cwd() / "LDA_2D_Scatter.html")


        lda_3dot = px.scatter_3d(mmm_lda, x = 'x', y = 'y', z = 'z', color = "Set", symbol = "Set", title = f"LDA of MMM 3D", hover_data = {"Fusion": True})
        lda_3dot.show()
        lda_3dot.write_html(pathlib.Path.cwd() / "LDA_3D_Scatter.html")




    def meta_stats(self, dataframe: pandas.DataFrame, name: str):

        dataframe = dataframe.unstack()
        # print(dataframe)
        # print(dataframe.mean())
        meta_data = pandas.Series(data = {"Mean": dataframe.mean(), "Median": dataframe.median(), "Mode": [entry for entry in dataframe.mode()], "Max": dataframe.max(), "Min": dataframe.min(), "SD": dataframe.std(), "Var": dataframe.var(), "Skew": dataframe.skew(), "Kurtosis": dataframe.kurtosis()}, name = name)
            
        return meta_data



    def pca_output(self):
        '''
        Because right now I don't know how to handle the mode (of if I even need to) I'm going to manually strip out the columns that I want. 
        '''
        if self.pca_output_on:


            mmm_frame = self.convert_2_R(self.mmm_frame)
            mtotal_fre = self.convert_2_R(self.Frequency.T)
            # print(robjects.r.rownames(mtotal_fre))

            pcaR = self.r_pca(mmm_frame)
            pcaF = self.r_pca(mtotal_fre)

            # This line does not work... don't know why.
            # self.pcrR.loadings = self.pcrR.rx2("rotation") * -1.0

            pca_sta_plot = pcaR.rx2("x")
            pca_fre_plot = pcaF.rx2("x")
            # print(self.pcr_plot)
            # print(type(self.pcr_plot))
            # sanity = self.pcr_plot.rx2("PC1")

            pca_sta_plot = robjects.vectors.DataFrame({'x': pca_sta_plot.rx(True, 1), 'y': pca_sta_plot.rx(True, 2), 'z': pca_sta_plot.rx(True, 3)})
            pca_fre_plot = robjects.vectors.DataFrame({'x': pca_fre_plot.rx(True, 1), 'y': pca_fre_plot.rx(True, 2), 'z': pca_fre_plot.rx(True, 3)})

            with localconverter(robjects.default_converter + pandas2ri.converter):
                pca_sta_plot: pandas.DataFrame = robjects.conversion.rpy2py(pca_sta_plot)
                pca_fre_plot: pandas.DataFrame = robjects.conversion.rpy2py(pca_fre_plot)
            
            pca_sta_plot["Set"] = pandas.Categorical(self.setnames)
            pca_fre_plot["Set"] = pandas.Categorical(["33P", "33Ran", "55P", "55Ran"])

            pca_sta_plotly = px.scatter_3d(pca_sta_plot, x = 'x', y = 'y', z = 'z', color = "Set", symbol = "Set", title = f"PCA of Statistics Data")
            pca_sta_plotly.show()
            pca_sta_plotly.write_html(pathlib.Path.cwd() / "PCA_3D_Scatter.html")

            pca_fre_plotly = px.scatter_3d(pca_fre_plot, x = 'x', y = 'y', z = 'z', color = "Set", symbol = "Set", title = f"PCA of Frequency Transposed")
            pca_fre_plotly.show()
            pca_fre_plotly.write_html(pathlib.Path.cwd() / "PCA_2D_Scatter.html")




    def r_pca(self, rataframe: robjects.vectors.DataFrame):
        '''
        '''

        sd_fix = '''
        function(rataframe) {
            rataframe <- rataframe[, apply(rataframe, 2, function(x) sd(x) > 0)]
            return(rataframe)
        }
        '''

        try:
            rataframe = robjects.r.prcomp(rataframe, scale = True, center = True)
        except RuntimeError:

            rsdfun = robjects.r(sd_fix)
            rataframe = rsdfun(rataframe)
            rataframe = robjects.r.prcomp(rataframe, scale = True, center = True)

        except rpy2.rinterface_lib.embedded.RRuntimeError:

            # print("### It actually caught the error! ###")
            rsdfun = robjects.r(sd_fix)
            rataframe = rsdfun(rataframe)
            rataframe = robjects.r.prcomp(rataframe, scale = True, center = True)

        except Exception as e:
            print(type(e))
            print("$$$ still moving to the other expcetions... $$$")

            rsdfun = robjects.r(sd_fix)
            rataframe = rsdfun(rataframe)
            rataframe = robjects.r.prcomp(rataframe, scale = True, center = True)

        return rataframe


    def random_comparison(self, agene, bgene):
        '''
        This finds the G Number of the a and b gene (note: a does not have to be the head nor does b have to be the tail)
        Then this finds a random sequence to comapre to.
        Finally this compares the G Number to the Random G Number via G' (setting the number to be 100 if RG = 0)
        '''
        g_num = self.g_number(agene, bgene)

        rgene = self.random_seq(len(agene))
        rg_num = self.g_number(agene, rgene)

        with numpy.errstate(divide='ignore'):
            g_prime = numpy.where(rg_num != 0., g_num / rg_num, 100)            

        return g_num, rg_num, g_prime



    def random_seq(self, length):
        seq = f''
        for _ in range(length):
            n = random.randint(1, 4)
            if n == 1:
                n = 'A'
            elif n == 2:
                n = 'C'
            elif n == 3:
                n = 'G'
            elif n == 4:
                n = 'T'
            seq = f'{seq}{n}'

        return seq

    
    def weights(self):
        '''
        Initializes the weights and another important matrix
        '''
        self.weight = numpy.array([(1 / (2**(i))) for i in range(1, self.length + 1)])

        # Generates a list of all possible Dissimilarity Numbers for this length
        gnums = [numpy.reshape(numpy.array(i), (1, self.length)) for i in itertools.product([0, 1], repeat = 1*self.length)]
        # Creates an emtpy dataframe for the running total of frequencies for each Dissimilarity Number
        zero_data = [0 for _ in range(0, 2**self.length)]

        self.Frequency = pandas.DataFrame(data = {colname: zero_data for colname in self.master_freq}, index = [numpy.sum(g * self.weight) for g in gnums])



    def write2excel(self):
        '''
        '''
        if self.write_files:
            if not self.outputpath.is_dir():
                self.outputpath.mkdir()
                # Can give a FileNotFoundError if parents don't exist: maybe add try except to create the path?
                # Can also give FileExistsError if the directory already is there, but I think this should be fine


            with pandas.ExcelWriter(self.outputpath / f"FFP_Matrix_Length-{self.length}.xlsx") as FFP_writer, \
                pandas.ExcelWriter(self.outputpath / f"TTP_Matrix_Length-{self.length}.xlsx") as TTP_writer, \
                pandas.ExcelWriter(self.outputpath / f"FFR_Matrix_Length-{self.length}.xlsx") as FFR_writer, \
                pandas.ExcelWriter(self.outputpath / f"TTR_Matrix_Length-{self.length}.xlsx") as TTR_writer, \
                pandas.ExcelWriter(self.outputpath / f"Similarity_Frequency_Length-{self.length}.xlsx") as Fre_writer, \
                pandas.ExcelWriter(self.outputpath / f"Statistics_Length-{self.length}.xlsx") as Sta_writer:

                self.Frequency.to_excel(Fre_writer, "Totals")

                for sheetname, _ in self.FFP_xlsx.items():
                    self.FFP_xlsx[sheetname].to_excel(FFP_writer, sheetname)
                    self.TTP_xlsx[sheetname].to_excel(TTP_writer, sheetname)
                    self.FFR_xlsx[f"{sheetname}_Random"].to_excel(FFR_writer, f"{sheetname}_Random")
                    self.TTR_xlsx[f"{sheetname}_Random"].to_excel(TTR_writer, f"{sheetname}_Random")
                    self.Fre_xlsx[sheetname].to_excel(Fre_writer, sheetname)
                
                for sheetname, _ in self.Stat_xlsx.items():
                    self.Stat_xlsx[sheetname].to_excel(Sta_writer, sheetname)


if __name__ in '__main__':
    ee = EE()