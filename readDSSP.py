import re
import sys

def sa2c(sa,aa):
    #maximum solvent accessiblity
    maxsa = {'A' : 106,
             'B' : 160,
             'C' : 135,
             'D' : 163,
             'E' : 194,
             'F' : 197,
             'G' : 84,
             'H' : 184,
             'I' : 169,
             'K' : 205,
             'L' : 164,
             'M' : 188,
             'N' : 157,
             'P' : 136,
             'Q' : 198,
             'R' : 248,
             'S' : 130,
             'T' : 142,
             'V' : 142,
             'W' : 227,
             'X' : 180,
             'Y' : 222,
             'Z' : 196}
    if re.compile(r"[a-z]").match(aa): return "F" # disulphide bridge
    if not maxsa(aa): return "-" # no amino acid
    rsa = sa / maxsa(aa)
    print("aa=%s sa=%5.1f  max_sa=%5.1f  rsa=%5.3f\n" %(aa, sa, maxsa(aa), rsa))
    if rsa <= 0.02: return "A"
    elif rsa <= 0.14: return "B"
    elif rsa <= 0.33: return "C"
    elif rsa <= 0.55: return "D"
    else: return "E"

def readDSSP(dsspfile, qrange):
    aa_dssp = ""
    sa_dssp = ""
    ss_dssp = ""
    
    open_dsspfile = open(dsspfile)
    dsspfile_lines = open_dsspfile.readlines()
    for line in dsspfile_lines:
        if re.compile(r"^\s*\#\s*RESIDUE\s+AA").match(line):
            break
        if re.compile(r"^.{5}(.{5})(.)(.)\s(.).\s(.).{18}(...)").match(line):
            l = line.split()
            thisres, icode, chain, aa, ss, sa = l[:6]
            range = qrange
            
            if aa == "!": continue
            
            thisres, chain, icode, sa = [int(x) for x in [thisres, chain, icode, sa]]
            
            contained = 0
            
            if qrange != '':
                while contained == 0 & range != "":
                    #syntax (A:56S-135S)
                    if re.sub(r"^(\S):(-?\d+)[A-Z]-(\d+)([A-Z])", "", range) and \
                        chain == l[0] and icode == l[3] and \
                            float(l[1]) <= thisres and thisres <= float(l[2]):
                                contained = 1
                    #syntax (R:56-135)
                    elif re.sub(r"^(\S):(-?\d+)[A-Z]?-(\d+)[A-Z]?", "", range) and \
                        chain == l[0] and l[1] <= thisres and thisres <= l[2]:
                            contained = 1
                    #syntax (56-135)
                    elif re.sub(r"^(-?\d+)[A-Z]-(\d+)([A-Z])", '', range) and \
                        chain == '' and icode == l[2] and l[0] <= thisres and thisres <= l[1]:
                            contained = 1
                    #syntax (56-135)
                    elif re.sub(r"^(-?\d+)[A-Z]?-(\d+)[A-Z]?", "",range) and \
                        chain == '' and l[0] <= thisres and thisres <= l[1]:
                            contained = 1
                    #syntax (A:) or (A:,2:)
                    elif re.sub(r"^(\S):", '', range) and chain == l[0]:
                        contained = 1
                    #syntax (-)
                    elif re.sub(r"^-$", "", range) and chain == '':
                        contained = 1

                    range = re.sub(r"^,", "", range)
                    
                if contained == 0:
                    continue
            aa_dssp += aa
            ss_dssp += ss
            sa_dssp += sa2c(sa,aa)
    open_dsspfile.close()
    return aa_dssp, ss_dssp, sa_dssp