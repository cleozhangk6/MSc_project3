import pandas as pd
import sys

# Create Dataframe for domain info
# Create dataframe
domains = []
starts = []
ends = []
colors = []

summary_file = sys.argv[1]
input_file = sys.argv[2]
with open(summary_file) as file:
# with open('P43238_summary.txt') as file:
    for line in file.readlines():
        cols = line.split()
        domains.append(cols[3]) # domain name
        starts.append(cols[7]) # start position
        ends.append(cols[9]) # end position

domains_af = pd.DataFrame({
        'Domain': domains,
        'Start': starts,
        'End': ends, 
        'Color': colors},
        index=None)



# # Create .pml commands to color domain
# with open('colbydom.pml','w') as f:
#     f.write('color gray, all')
#     i = 2
#     for index, row in domains_af.iterrows():
#         drange = str(row['Start']) + '-' + str(row['End'])
#         dname = row['Domain']
#         f.write('; create ' + dname + ',resi ' + drange + 
#         '; color ' + str(i) + ',' + dname)
#         i += 1 if i!=6 else 2  # skip as 6 and 7 same color

#     f.write('\n zoom' +
#             '\n bg_color white' + 
#             '\n viewport 800, 800' +
#             '\n set depth_cue,1' +
#             '\n set antialias, 1' +
#             '\n set cartoon_fancy_sheets, 1' +
#             '\n set two_sided_lighting, on' +
#             '\n set cartoon_loop_quality=100' +
#             '\n set ray_opaque_background, on' +
#             '\n set ray_shadows,off' + 
#             '\n ray')

pymol_colors = [('white', 0), ('black', 1), ('blue', 2), ('green', 3), ('red', 4), ('cyan', 5), ('yellow', 6), ('dash', 7), ('magenta', 8), ('salmon', 9), ('lime', 10), ('slate', 11), ('hotpink', 12), ('orange', 13), ('chartreuse', 14), ('limegreen', 15), ('purpleblue', 16), ('marine', 17), ('olive', 18), ('purple', 19), ('teal', 20), ('ruby', 21), ('forest', 22), ('deepblue', 23), ('grey', 24), ('gray', 25), ('carbon', 26), ('nitrogen', 27), ('oxygen', 28), ('hydrogen', 29), ('brightorange', 30), ('sulfur', 31), ('tv_red', 32), ('tv_green', 33), ('tv_blue', 34), ('tv_yellow', 35), ('yelloworange', 36), ('tv_orange', 37), ('pink', 48), ('firebrick', 49), ('chocolate', 50), ('brown', 51), ('wheat', 52), ('violet', 53)]

# Create .spt file for JSmol
with open('../myJmolscript.spt','w') as f:
    f.write('load uploads/' + input_file + '; cartoon ONLY')
    i = 2
    for index, row in domains_af.iterrows():
        drange = str(row['Start']) + '-' + str(row['End'])
        dname = row['Domain']
        row['Color'] = pymol_colors[i][0]
        f.write(
        '; select ' + drange + 
        '; color ' + pymol_colors[i][0])
        i += 1 if i!=6 else 2  # skip as 6 and 7 same color


# domains_af.to_csv('domain_info.txt', header=None, index=None, sep=' ', mode='a')
domains_af.to_json('domains_info.json', orient="records")