import pandas as pd
import sys

# Create Dataframe for domain info
# Create dataframe
domains = []
starts = []
ends = []
filename = sys.argv[1]
with open(filename) as file:
# with open('P43238_summary.txt') as file:
    for line in file.readlines():
        cols = line.split()
        domains.append(cols[3]) # domain name
        starts.append(cols[7]) # start position
        ends.append(cols[9]) # end position

domains_af = pd.DataFrame({
        'Domain': domains,
        'Start': starts,
        'End': ends})

print(domains_af)

# Create .pml commands to color domain
with open('colbydom.pml','w') as f:
    f.write('color gray, all')
    i = 2
    for index, row in domains_af.iterrows():
        drange = str(row['Start']) + '-' + str(row['End'])
        dname = row['Domain']
        f.write('; create ' + dname + ',resi ' + drange + 
        '; color ' + str(i) + ',' + dname)
        i += 1 if i!=6 else 2  # skip as 6 and 7 same color

    f.write('\n zoom' +
            '\n bg_color white' + 
            '\n viewport 1600, 1600' +
            '\n set depth_cue,1' +
            '\n set antialias, 1' +
            '\n set cartoon_fancy_sheets, 1' +
            '\n set two_sided_lighting, on' +
            '\n set cartoon_loop_quality=100' +
            '\n set ray_opaque_background, on' +
            '\n set ray_shadows,off' +
            '\n ray')