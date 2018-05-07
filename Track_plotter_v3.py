import ROOT, math, sys, os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
####################################################################################################
def read_file(filename):
    with open(filename) as file:
        return file.readlines()

def make_dir(dir_name):
	if(os.path.exists(dir_name) != True):
		os.system('mkdir ' + dir_name)

#Parameters
nor_size = 5
col_size = 30
fa = 16
alpha = 5
phi = 25
evr = 90
theta = 148


study_dir = sys.argv[2][:-5] + "_Plots/"
#p = ROOT.TFile.Open("results_alpha_%d_ConeAng_%d_evr_0_%d_theta_%d.root"%(alpha,phi,evr,theta))
p = ROOT.TFile.Open(sys.argv[2])

make_dir(study_dir)

#What to plot
#option = 'c'
option = sys.argv[1]

if option == 'a':
	bool_michels = True
	bool_selected = True
elif option == 'b':
	bool_michels = False
	bool_selected = True
elif option == 'c':
	bool_michels = True
	bool_selected = False
else:
	sys.exit()

#a_michels_sel
if bool_michels and bool_selected:
	dir_name = study_dir + "a_michels_selected/"
	michel_find = 1
	sel_find = 1
elif bool_michels != True and bool_selected:
	dir_name = study_dir + "b_not_michels_selected/"
	michel_find = 0
	sel_find = 1
elif bool_michels and bool_selected != True:
	dir_name = study_dir + "c_michels_not_selcted/"
	michel_find = 1
	sel_find = 0
else:
	sys.exit()

#Read in CSV Michel file
vertex_rn = []
vertex_ev = []
vertex_cl = []
vertex_x = []
vertex_y = []
vertex_z = []

csv_michel = "Michel_candidates_vertex_v8.csv"
csv_count = 0
with open(csv_michel) as csvfile:
	reader = csv.reader(csvfile)
	for row in reader:
		if csv_count == 0:
			csv_count += 1
			continue
		#print(row[2])
		vertex_rn.append(float(row[0]))
		vertex_ev.append(float(row[1]))
		vertex_cl.append(float(row[2]))
		vertex_x.append(float(row[3]))
		vertex_y.append(float(row[4]))
		vertex_z.append(float(row[5]))

#make needed directories
make_dir(dir_name)

#read in root file
michel_run = []
michel_eve = []
michel_trk = []


lalpha_x = []
lalpha_y = []
lalpha_z = []
ford_x = []
ford_y = []
ford_z = []
dist = []
for point in p.nt_brk_pts:
	michel_run.append(point.run_num)
	michel_eve.append(point.ev_num)
	michel_trk.append(point.cluster_id)
	lalpha_x.append(point.lalpha_x)
	lalpha_y.append(point.lalpha_y)
	lalpha_z.append(point.lalpha_z)
	ford_x.append(point.ford_x)
	ford_y.append(point.ford_y)
	ford_z.append(point.ford_z)
	dist.append(point.dist)

pca_run = []
pca_eve = []
pca_trk = []
pca_x = []
pca_y = []
pca_z = []

for r in p.nt_pca_points:
	pca_run.append(r.run_num)
	pca_eve.append(r.ev_num)
	pca_trk.append(r.cluster_id)
	pca_x.append(r.pca_x)
	pca_y.append(r.pca_y)
	pca_z.append(r.pca_z)

vec_run = []
vec_eve = []
vec_trk = []
vec_s_x = []
vec_s_y = []
vec_s_z = []
vec_pca_x =[]
vec_pca_y =[]
vec_pca_z =[]
vec_np_x = []
vec_np_y = []
vec_np_z = []
vec_angle = []

for v in p.nt_vectors:
	vec_run.append(v.run_num)
	vec_eve.append(v.ev_num)
	vec_trk.append(v.cluster_id)
	vec_s_x.append(v.pca_x_s)
	vec_s_y.append(v.pca_y_s)
	vec_s_z.append(v.pca_z_s)
	vec_pca_x.append(v.pca_x_e)
	vec_pca_y.append(v.pca_y_e)
	vec_pca_z.append(v.pca_z_e)
	vec_np_x.append(v.lalpha_x_e)
	vec_np_y.append(v.lalpha_y_e)
	vec_np_z.append(v.lalpha_z_e)
	vec_angle.append(v.angle)

count = 0
pts_run = []
pts_eve = []
pts_clu = []
pts_x = []
pts_y = []
pts_z = []
pts_mch = []
pts_sel = []
pts_ord =[]

for m in p.nt_trk_pts:
	if(m.michel == michel_find and m.selected == sel_find):
		old_suma = str(m.run_num) + str(m.ev_num) + str(m.cluster_id)
		break	

all_clusters = []
n_cluster = []

end_trk_points = p.nt_trk_pts.GetEntries()

count = 0
for pts in p.nt_trk_pts:
	count +=1
	suma = str(pts.run_num) + str(pts.ev_num) + str(pts.cluster_id)
	if(pts.michel == michel_find and pts.selected == sel_find):
		if suma == old_suma:
			a = []
			"""
			pts_run.append(pts.run_num)
			pts_eve.append(pts.ev_num)
			pts_clu.append(pts.cluster_id)
			pts_x.append(pts.x)
			pts_y.append(pts.y)
			pts_z.append(pts.z)
			pts_mch.append(pts.michel)
			pts_sel.append(pts.selected)
			pts_ord.append(pts.point_ord)
			"""
			a.append(pts.run_num)
			a.append(pts.ev_num)
			a.append(pts.cluster_id)
			a.append(pts.x)
			a.append(pts.y)
			a.append(pts.z)
			a.append(pts.michel)
			a.append(pts.selected)
			a.append(pts.point_ord)
			a.append(pts.vertex)
			n_cluster.append(a)
			old_suma = suma
		else:
			all_clusters.append(n_cluster)
			a = []
			n_cluster = []
			a.append(pts.run_num)
			a.append(pts.ev_num)
			a.append(pts.cluster_id)
			a.append(pts.x)
			a.append(pts.y)
			a.append(pts.z)
			a.append(pts.michel)
			a.append(pts.selected)
			a.append(pts.point_ord)
			a.append(pts.vertex)
			n_cluster.append(a)
			old_suma = suma
		if(count == end_trk_points-1):
			all_clusters.append(n_cluster)			
print(len(all_clusters))

break_counter = 0
for cluster in all_clusters:
#	if break_counter == 2:
#		break
	break_counter += 1
	#print(break_counter)
	#Finding Vectors
	for v in xrange(len(vec_pca_x)):
		if(vec_run[v] == cluster[0][0] and vec_eve[v] == cluster[0][1] and vec_trk[v] == cluster[0][2]):
			soa_xy = np.array([[vec_s_x[v],vec_s_y[v],fa*vec_pca_x[v],fa*vec_pca_y[v]], [vec_s_x[v],vec_s_y[v],fa*vec_np_x[v],fa*vec_np_y[v]]])
			soa_zy = np.array([[vec_s_z[v],vec_s_y[v],fa*vec_pca_z[v],fa*vec_pca_y[v]], [vec_s_z[v],vec_s_y[v],fa*vec_np_z[v],fa*vec_np_y[v]]])
			soa_zx = np.array([[vec_s_z[v],vec_s_x[v],fa*vec_pca_z[v],fa*vec_pca_x[v]], [vec_s_z[v],vec_s_x[v],fa*vec_np_z[v],fa*vec_np_x[v]]])
			angle = vec_angle[v]
			angle_xy = math.acos(vec_pca_x[v]*vec_np_x[v] + vec_pca_y[v]*vec_np_y[v]) * 180./3.14159
			angle_zy = math.acos(vec_pca_z[v]*vec_np_z[v] + vec_pca_y[v]*vec_np_y[v]) * 180./3.14159
			angle_zx = math.acos(vec_pca_z[v]*vec_np_z[v] + vec_pca_x[v]*vec_np_x[v]) * 180./3.14159
			break
		else:
			soa_xy = np.array([[0,0,0,0], [0,0,0,0]])
			soa_yz = np.array([[0,0,0,0], [0,0,0,0]])
			soa_zx = np.array([[0,0,0,0], [0,0,0,0]])
			angle = 0.
	x_ord = []
	y_ord = []
	z_ord = []
	x_un = []
	y_un = []
	z_un = []
	col_ord = []
	size_ord = []
	pca_count_ord = 0
	found_point_r_ord = False
	found_point_y_ord = False
	found_pca_ord = False
	found_vertex_g = False
	print("###############")
	print"%d, %d, %d"%(cluster[0][0],cluster[0][1],cluster[0][2])
	print(len(cluster))
	for point in cluster:
		found_point_r_ord = False
		found_point_y_ord = False
		found_pca_ord = False
		found_vertex_g = False
		found_hyp_vertex = False
		#print"%d, %d. %d, %f, %f, %f, %d, %d. %d"%(point[0],point[1],point[2],point[3],point[4],point[5], point[6], point[7], point[8])
		if point[8] == 1:
			x_ord.append(point[3])
			y_ord.append(point[4])
			z_ord.append(point[5])
			#print"%d, %d. %d, %f, %f, %f"%(point[0],point[1],point[2],point[3],point[4],point[5]
			for cut_n in xrange(len(vertex_x)):
				if(vertex_rn[cut_n] == point[0] and vertex_ev[cut_n] == point[1] and vertex_cl[cut_n] == point[2] and vertex_x[cut_n] + .001 > point[3] and vertex_x[cut_n] - .001 < point[3] and vertex_y[cut_n] + .001 > point[4] and vertex_y[cut_n] - .001 < point[4] and vertex_z[cut_n] + .001 > point[5] and vertex_z[cut_n] - .001 < point[5]):
					print"%f, %f, %f"%(point[3],point[4],point[5])
					vert_x_point_ord =  point[3]
					vert_y_point_ord =  point[4]
					vert_z_point_ord =  point[5]
					col_ord.append('g')
					found_vertex_g = True
					size_ord.append(col_size+8)
					break
			if(found_vertex_g):
				continue
			if bool_selected != True:
				for cut_n in xrange(len(ford_x)):
					if(ford_x[cut_n] + .001 > point[3] and ford_x[cut_n] - .001 < point[3] and ford_y[cut_n] + .001 > point[4] and ford_y[cut_n] - .001 < point[4] and ford_z[cut_n] + .001 > point[5] and ford_z[cut_n] - .001 < point[5]):
						ford_x_point_ord =  point[3]
						ford_y_point_ord =  point[4]
						ford_z_point_ord =  point[5]
						col_ord.append('r')
						found_point_r_ord = True
						distance_ord = dist[cut_n]
						size_ord.append(col_size+5)
						break
				for i in xrange(len(pca_x)):
					if(pca_run[i] == point[0] and pca_eve[i] == point[1] and pca_trk[i] == point[2] and pca_x[i] + .001 > point[3] and pca_x[i] - .001 < point[3] and pca_y[i] + .001 > point[4] and pca_y[i] - .001 < point[4] and pca_z[i] + .001 > point[5] and pca_z[i] - .001 < point[5]):
						#print(pca_count_ord)
						found_pca_ord = True
						if pca_count_ord == 9:
							col_ord.append('b')
							size_ord.append(col_size)
							#print("Last BLUE")
							break
						else:
							col_ord.append('b')
							size_ord.append(col_size)
							#print("BLUE")
						pca_count_ord += 1	
			if point[9] == 1:
				col_ord.append('c')
				found_hyp_vertex = True
				size_ord.append(col_size+10)
			
			
			if(found_point_r_ord != True and found_point_y_ord != True and found_pca_ord != True and found_vertex_g != True and found_hyp_vertex != True):
				col_ord.append('k')
				size_ord.append(nor_size)	
	x_un = []
	y_un = []
	z_un = []
	col_un = []
	size_un = []
	
	for point in cluster:
		found_vertex_g_un = False
		found_point_y_un = False
		#print"%d, %d. %d, %f, %f, %f, %d, %d. %d"%(point[0],point[1],point[2],point[3],point[4],point[5], point[6], point[7], point[8])
		if point[8] == 0:
			x_un.append(point[3])
			y_un.append(point[4])
			z_un.append(point[5])
			#print"%d, %d. %d, %f, %f, %f"%(point[0],point[1],point[2],point[3],point[4],point[5]
			if bool_selected != True:
				for cut_n in xrange(len(ford_x)):
					if(lalpha_x[cut_n] + .001 > point[3] and lalpha_x[cut_n] - .001 < point[3] and lalpha_y[cut_n] + .001 > point[4] and lalpha_y[cut_n] - .001 < point[4] and lalpha_z[cut_n] + .001 > point[5] and lalpha_z[cut_n] - .001 < point[5]):
						lalpha_x_point_un =  point[3]
						lalpha_y_point_un =  point[4]
						lalpha_z_point_un =  point[5]
						col_un.append('y')
						found_point_y_un = True
						size_un.append(col_size+5)
						break
			for cut_n in xrange(len(vertex_x)):
				if(vertex_rn[cut_n] == point[0] and vertex_ev[cut_n] == point[1] and vertex_cl[cut_n] == point[2] and vertex_x[cut_n] + .001 > point[3] and vertex_x[cut_n] - .001 < point[3] and vertex_y[cut_n] + .001 > point[4] and vertex_y[cut_n] - .001 < point[4] and vertex_z[cut_n] + .001 > point[5] and vertex_z[cut_n] - .001 < point[5]):
					vert_x_point_un =  point[3]
					vert_y_point_un =  point[4]
					vert_z_point_un =  point[5]
					col_un.append('g')
					found_vertex_g_un = True
					size_un.append(col_size+5)
					break
			if(found_point_y_un != True and found_vertex_g_un != True):
				col_un.append('m')
				size_un.append(nor_size)
	print"Ordered tracks = %d"%(len(col_ord))
	print"Unordered points = %d"%(len(col_un))
	#print(len(x))
	XYZplots = plt.figure()
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	if bool_selected != True:
		X, Y, U, V = zip(*soa_xy)
		XYplot.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1, color = ['c','y'])
	XYplot.scatter(x_ord ,y_ord, s = size_ord, c = col_ord,linewidth=0)
	XYplot.plot(x_ord ,y_ord)
	XYplot.scatter(x_un ,y_un, s = size_un, c = col_un,linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)
	if bool_selected != True:
		X, Y, U, V = zip(*soa_zy)
		ZYplot.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1, color = ['c','y'])
	ZYplot.scatter(z_ord, y_ord, s = size_ord, c = col_ord,linewidth=0)
	ZYplot.plot(z_ord,y_ord)
	ZYplot.scatter(z_un ,y_un, s = size_un, c = col_un,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	if bool_selected != True:
		X, Y, U, V = zip(*soa_zx)
		XZplot.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1, color = ['c','y'])
	XZplot.scatter(z_ord, x_ord, s = size_ord, c = col_ord,linewidth=0)
	XZplot.plot(z_ord, x_ord)
	XZplot.scatter(z_un ,x_un, s = size_un, c = col_un,linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	Text_box = XYZplots.add_subplot(224)
	Text_box.axis([0, 10, 0, 10])
	Text_box.set_xticks([])
	Text_box.set_yticks([])
	Track_info = 'Run: %d, Event: %d, Cluster: %d'%(point[0],point[1],point[2])
	File_name = 'Run_%d_Event_%d_Cluster_%d'%(point[0],point[1],point[2])
	#Info = 'Alpha = %d; Distance = %f'%(alpha, distance)
	#ford_info = 'Final point in ordered cluster:'\
	#			'(%f, %f, %f)'%(ford_x_point,ford_y_point,ford_z_point)
	#lalpha_info = 'Next point tested in cluster: '\
	#			'(%f, %f, %f)'%(lalpha_x_point,lalpha_y_point,lalpha_z_point)
	angle_info = 'Angle between vectors = %f'%(angle)
	#angle_info2 = 'XY angle = %d; ZY angle = %d; ZX angle = %d'%(angle_xy,angle_zy,angle_zx)
	Text_box.text(9, 9, Track_info, family='arial', ha='right', fontsize=13)
	#Text_box.text(9, 7, Info, family='arial', ha='right', fontsize=13)
	#Text_box.text(9.5, 5, ford_info, family='arial', ha='right', fontsize=8)
	#Text_box.text(9.5, 4.5, lalpha_info, family='arial', ha='right', fontsize=8)
	#if bool_selected != True:
	#	Text_box.text(9,8,angle_info,family='arial', ha='right', fontsize=13)
	#Text_box.text(8, 2, angle_info2,family='arial', ha='right', fontsize=8)
	XYZplots.set_size_inches(10,10)
	XYZplots.savefig(dir_name + File_name + ".pdf")
	XYZplots.clf()


sys.exit()