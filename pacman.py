#!/dls_sw/apps/python/anaconda/1.7.0/64/bin/python
import os, sys, time
import numpy as np
from Chip_StartUp_v5 import get_shot_order
from Chip_StartUp_v5 import get_alphanumeric
from Chip_StartUp_v5 import get_xy
from matplotlib import pyplot as plt

def hits_scrape(fid, column_choice, shot_order_addr_list, bound=False):
    f = open(fid)
    listy = f.readlines()[1].split('|')
    f.close()
    print 'COLUMN CHOICES'
    for i, item in enumerate(listy[2:10]):
        print 'Column', i+1, item

    hits_dict = {}
    f = open(fid)
    for i, line in enumerate(f.readlines()[3:]):
        if line.startswith('-'):
            break
        line.rstrip('\n')
        listy = line.split('|')
        entry = [x.lstrip().rstrip() for x in listy[1:]]
        if not entry[column_choice]:
            continue 
        addr = shot_order_addr_list[i]
        if bound:
            hits_dict[addr] = i 
        else:
            hits_dict[addr] = entry[column_choice] 
    f.close()
    return hits_dict 

def det_dist_file_scrape(fid, shot_order_addr_list):
    print shot_order_addr_list[:10]  
    
    det_dict = {} 
    f = open(fid)
    for line in f.readlines():
        entry = line.split()
        img_num = int(entry[0])
        val = float(entry[1])
        if val == 0:
            val = None
        det_dict[img_num] = val
    
    hits_dict = {}
    for i in range(len(shot_order_addr_list)):
        addr = shot_order_addr_list[i]
        if i in det_dict.keys():
            val = det_dict[i]
        else:
            val = None
        hits_dict[addr] = val
    return hits_dict

def make_plot_arrays(hits_addr_dict, chip_type):
    x_list, y_list, z_list = [],[],[]
    for addr, val in hits_addr_dict.items():
       x, y = get_xy(addr, chip_type)
       x_list.append(x)
       y_list.append(y)
       z_list.append(val)

    X = np.array(x_list)
    Y = np.array(y_list)
    Z = np.array(z_list)

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()
    return x, y, z 

def run_fromdir_method(*args):
    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
	if k.startswith('dir'):
	    path = v
            print '\n\n\nThis is the path:', path
        elif k.startswith('wild'):
            wild = str(v)
            print 'This is the wild:', path

    addr_list = get_shot_order(chip_type='1')
    hits_dict = {}
    if path:	
        if os.path.isdir(path):
            fid_list = os.listdir(path)
            int_list = [fid for fid in fid_list if fid.startswith('int-0')]
	    print int_list[:10]
            shotnum_list = [int(fid.split('_')[1].rstrip('.pickle')) for fid in int_list]
            print 'Number of Integrated:', len(shotnum_list)
	    for x in  shotnum_list[:20]:
                print x,

	    for i, addr in enumerate(addr_list): 
		if i in shotnum_list:
		    hits_dict[addr] = 1
		else:
		    hits_dict[addr] = 0
            x, y, z = make_plot_arrays(hits_dict, chip_type='1') 
            return x, y, z 
        else:
            print 'Not a Valid Path'
            return 0
    else:
        raise SyntaxError('Unknown directory\n\n\n')
        return 0
    
def run_fromfile_method(*args):
    fid = None
    col = 3
    bind_type = None
    chip_type = '1'
    block_list = []
    for arg in args:
        k = arg.split('=')[0]
        v = arg.split('=')[1]
        if k.startswith('file'):
            fid = v
        elif k.startswith('col'):
            col = int(v)
        elif k.startswith('chiptype'):
            chip_type = v
        elif k.startswith('bind'):
            bind_type = v
        elif k.startswith('blocks'):
            try:
                block_list = v.split(',')
                for block in block_list: print block
            except SyntaxError('Expected list of block separated by commas: A1,A2,A4'):
                return 0

    if fid:
        if bind_type:
            if bind_type.startswith('alpha'):
                addr_list = get_alphanumeric(chip_type)
            elif bind_type.startswith('shot'):
                addr_list = get_shot_order(chip_type)
            bound = True
        else:
            addr_list = get_shot_order(chip_type)
            bound = False

        print 'Length of address list:%s\n' %len(addr_list)
        new_addr_list = []
        if len(block_list) > 0:
            for addr in addr_list:
                block = addr.split('_')[0]
                if block in block_list:
                    new_addr_list.append(addr)
            addr_list = new_addr_list
            print 'After:            ', len(addr_list)

        if fid.endswith('out'):
            hits_addr_dict = hits_scrape(fid, col, addr_list, bound)           
            x, y, z = make_plot_arrays(hits_addr_dict, chip_type) 
            return x, y, z
        elif fid.endswith('.txt'):
            hits_addr_dict = det_dist_file_scrape(fid, addr_list)           
            x, y, z = make_plot_arrays(hits_addr_dict, chip_type) 
            return x, y, z
        else:
            raise SyntaxError('Expected file that ends with .out\n\n\n')
            return 0
            
    else:
        raise SyntaxError('Expected file that ends with .spots or .out')
        print 'Insert raw intput code here'
        return 0

def plot(x, y, z, *args):
    alfa = 1
    col = 3
    mrksz = 8
    pixels = 200
    cmap_choice = 'terrain'
    xlim_min, xlim_max = (-1,25)
    ylim_min, ylim_max = (-1,25)
    zlim_min, zlim_max = (None,None)

    for arg in args:
        k = arg.split('=')[0]
        v = arg.split('=')[1]
        if k.startswith('file'):
            fid = v
        if k.startswith('dir'):
            fid = v 
        elif k.startswith('col'):
            col = str(v)
        elif k.startswith('ms'):
            mrksz = int(v)
        elif k.startswith('alpha'):
            alfa = float(v)
        elif k.startswith('cmap'):
            cmap_choice = str(v)
        elif k.startswith('dpi'):
            pixels = int(v)
        elif k.startswith('xlim'):
            xlim_min, xlim_max = v.split(',')
        elif k.startswith('ylim'):
            ylim_min, ylim_max = v.split(',')
        elif k.startswith('zlim'):
            zlim_min, zlim_max = v.split(',')
            zlim_min = float(zlim_min)
            zlim_max = float(zlim_max)

    raw_fid = fid.split('/')[-1]
    print '\nraw_fid', raw_fid
    if raw_fid == '':
        raw_fid = fid.split('/')[-2]


    out_fid = '%s_%s_%smap.png' %(raw_fid.split('.')[0], col, time.strftime("%Y%m%d_%H%M%S"))
    fig = plt.figure(figsize=(10,10), facecolor='0.75', edgecolor='w')
    ax1 = fig.add_subplot(111, aspect=1, axisbg='0.5')
    fig.subplots_adjust(left=0.03,right=0.97,bottom=0.03,top=0.97)

    #plt.axis('off')
    ax1.set_xlim(float(xlim_min), float(xlim_max))
    ax1.set_ylim(float(ylim_min), float(ylim_max))
    
    if not (zlim_min, zlim_max) == (None, None):
        a_list = []
        for entry in z.tolist():
            if float(entry) > zlim_max:
                a_list.append(zlim_max)
            elif float(entry) < zlim_min:
                a_list.append(zlim_min)
            else:
                a_list.append(entry)
        z = np.array(a_list)

    plt.scatter(x, y, c=z, s=mrksz, marker='s', alpha=alfa, cmap=cmap_choice)

    ax1.invert_yaxis()
    cb = plt.colorbar()
    path = '/dls_sw/i24/scripts/pacman-pictures/'
    plt.savefig(path+out_fid, dpi=pixels, facecolor='0.15', bbox_inches='tight', pad_inches=0, transparent=True)
    print '\n\nPLOTTING KEYWORDS\n', 30*'-'
    print 'col=%s'     %col
    print 'ms=%s'      %mrksz
    print 'cmap=%s'    %cmap_choice
    print 'dpi=%s'     %pixels
    print 'xlim=%s,%s' %(xlim_min, xlim_max)
    print 'ylim=%s,%s' %(ylim_min, ylim_max)
    print 'zlim=%s,%s' %(zlim_min, zlim_max)
    print 'alpha=%s'   %alfa
    print 30*'-', '\n\n'
    print 'Imaged saved as: %s%s' %(path, out_fid)
    plt.show()

def run_dosecolumns(*args):
    path = '/dls/i24/data/2017/nt14493-54/processing/find_spots/hector_pacman/'
    chip_name = 'hector'
    num_of_doses = 5
    list_of_lists = []
    for x in range(num_of_doses):
        list_of_lists.append([])

    if os.path.isdir(path):
        fid_list = sorted(os.listdir(path))
        list_of_files = [fid for fid in fid_list if fid.startswith(chip_name)]
        for fid in list_of_files:
           print fid
           f = open(path+fid, 'r')
           for i, line in enumerate(f.readlines()[3:]):
               #print i, i % num_of_doses
               if line.startswith('-'):
                   break
               else:
                   list_of_lists[i % num_of_doses].append(line)
           f.close()

    f = open(path+fid, 'r')
    header = f.readlines()[:3]
    f.close()

    for i in range(num_of_doses):
        new_fid = chip_name + '.dose%s' %(i+1)
        print 'Making:', new_fid
        g = open(new_fid, 'a') 
        for line in header:
	    print line.rstrip()
            g.write(line)
        for line in list_of_lists[i]:
            g.write(line)
        g.close()
    print 
    
def main(*args):
    allowed_keyword_list = ['file','dir','binding','column','wild','chiptype', 'blocks',\
     'ms', 'cmap', 'dpi', 'xlim', 'ylim', 'zlim', 'alpha']

    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
        print "Keyword argument: %s=%s" % (k, v)
        there = [True for key in allowed_keyword_list if k in key]
        if True not in there:
            print 'Unknown arg in args:', arg
            print 'Allowed Keywords:', allowed_keyword_list
            print 'Exiting'
            return 0

    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
        if 'file' in k:
            method = 'fromfile'
        elif 'dir' in k:
            method = 'fromdir'
        elif 'dosecolumns' in k:
            method = 'dosecolumns'
    print "\n\nCalled PACMAN method = %s\n" % method
    print '\n\nDATA KEYWORDS\n', 30*'-'
    print 'file=filename.out'
    #print 'dir=/where/the/processed/data/is/kept/ ....... This for directory method NOT the file directory'
    #print 'dosecolumns=1'
    print 'blocks=A1,A2,A3,A4 ........................... Blocks MUST be in block order'
    print '                                               Total number of images in file must not exceed 200*num_of_blocks'
    print 30*'-', '\n'

    if method == 'fromfile':
        x, y, z = run_fromfile_method(*args)
    elif method == 'fromdir':
        x, y, z = run_fromdir_method(*args)
    elif method == 'dosecolumns':
        something = run_dosecolumns(*args)
        print 'Exiting'
        return 0 
    else:
        raise SyntaxError('Unknown method')
    
    plot(x, y, z, *args)
    print 'EOP'

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print '\n\t\t\t\tPACMAN'
        print './pacman.py file=/path/to/your/Filename.out'
        #print './pacman.py dir=YourDirectoryHere wild=yourwildcard\n'
        print "\t\t\t\tEXAMPLE\n./pacman.py file=egbert_spots.out col=3 ms=8 cmap='terrain' xlim=0,25 ylim=0,25\n"
        print '\t\t\t\tDOCUMENT\nsee https://docs.google.com/document/d/12564wj9aLwx7OlwqWHlEmI4lT_1qErsuvedVgjUE8F0/edit\n\n'
    else:
        main(*sys.argv[1:])
