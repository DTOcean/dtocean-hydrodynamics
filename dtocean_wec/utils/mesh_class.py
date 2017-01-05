import numpy as np
from numpy import array, zeros
import os

class Mesh():
    def __init__(self, path, fn):
        self.fn = os.path.join(path, fn)
        self.path = path
        if not os.path.isfile(self.fn):
            raise IOError("The filename or path are incorrect")
        mesh_type = self.fn[-3:]
        
        if mesh_type.lower() == "gdf":
            v, p = self.readWAMIT(self.fn)
            p -= 1
        elif mesh_type.lower() == "dat":
            # try to see if the specified filename is Nemoh format compliant
            v, p = self.readNemoh(self.fn)
        else:
            raise IOError('Mesh type not understood')
        
        self.v = v
        self.p = p
        n, c = self.__evalnorm()
        self.n = n
        self.c = c
        
    def __evalnorm(self):
        pan = self.p+1
        ver = self.v
        v12 = ver[pan[:,1]-1]-ver[pan[:,0]-1]
        v23 = ver[pan[:,2]-1]-ver[pan[:,1]-1]
        v34 = ver[pan[:,3]-1]-ver[pan[:,2]-1]
        v41 = ver[pan[:,0]-1]-ver[pan[:,3]-1]
        
        n1 = np.cross(v12, -v41)
        n2 = np.cross(v23, -v12)
        n3 = np.cross(v34, -v23)
        n4 = np.cross(v41, -v34)
        
            
        # only point 1 and 4 are repeated for the case of triangle
        mask_12 = np.linalg.norm(v12, axis=1) < 1e-6
        weight_12 = np.ones((n1.shape[0],1))
        weight_12[mask_12] = 0
        mask_23 = np.linalg.norm(v23, axis=1) < 1e-6
        weight_23 = np.ones((n1.shape[0],1))
        weight_23[mask_23] = 0
        mask_34 = np.linalg.norm(v34, axis=1) < 1e-6
        weight_34 = np.ones((n1.shape[0],1))
        weight_34[mask_34] = 0
        mask_41 = np.linalg.norm(v41, axis=1) < 1e-6
        weight_41 = np.ones((n1.shape[0],1))
        weight_41[mask_41] = 0
        
        
        mean_w = 4.0*np.ones((n1.shape[0],1))
        mean_w[mask_41] = 3.0
        mean_w[mask_12] = 3.0
        mean_w[mask_23] = 3.0
        mean_w[mask_34] = 3.0
        
        centroid = (ver[pan[:,0]-1]*weight_12+
                    ver[pan[:,1]-1]*weight_23+
                    ver[pan[:,2]-1]*weight_34+
                    ver[pan[:,3]-1]*weight_41)/mean_w
                    
            
        # keep the area as weight not the mean value,
        # the triangle (repeated) vertex are not adding 
        # value to the mean
#        n = (n1 + n2 + n3 + n4)/mean_w

        if mean_w[-1] == 3.0:
            n = n1/2.0
        else:
            n = (n1 + n3)/2.0
#        n /= np.linalg.norm(n)
#        print(np.linalg.norm(n))  
        return n, centroid
        
    def readNemoh(self, f_n):
        with open(f_n,'r') as mesh_f: 
            first_line = array(mesh_f.readline().split(),dtype = int)
            if len(first_line) == 1:
                nV = first_line
                nP = array(mesh_f.readline().split(),dtype = int)
                Vertex = zeros((nV,3))
                for vertex in range(nV):
                    Vertex[vertex,:] = array(mesh_f.readline().split(),dtype=float)
                Connectivity = zeros((nP,4))
                for panel in range(nP):
                    Connectivity[panel,:] = array(mesh_f.readline().split(),dtype=float)
            elif len(first_line) == 2:
                self.xsim = first_line[1]
                Vertex = []
                Connectivity = []
                pass_to_connectivity = False
                for line in mesh_f:
                    temp = array(line.split(),dtype = float)
                    if not int(temp[0]) == 0 and not pass_to_connectivity:
                        Vertex.append(temp[1:].tolist())
                    else:
                        pass_to_connectivity = True
                        Connectivity.append(temp.tolist())
                    
                Vertex = array(Vertex)
                Connectivity.pop(0)
                Connectivity.pop(-1)
                Connectivity = array(Connectivity,dtype=int)
            else:
                raise IOError('Mesh type not understood')
                
        return Vertex, Connectivity
        
    def readWAMIT(self, f_n):
        
        with open(f_n,'r') as mesh_f: 
                burnheader= mesh_f.readline(); del(burnheader)
                burnheader= mesh_f.readline(); del(burnheader)
                xsim= array(mesh_f.readline().split()[0],dtype=int)
                nP = array(mesh_f.readline().partition('!')[0],dtype=int)
                Vertex = zeros((nP*4,3))
                for vertex in range(nP*4):
                    Vertex[vertex,:] = array(mesh_f.readline().partition('!')[0].split(),dtype=float)
                Connectivity = zeros((nP,4),dtype=int)
                vertex = 1
                for panel in range(nP):
                    Connectivity[panel,:] = array([vertex,vertex+1,vertex+2,vertex+3])
                    vertex +=4
                    
        return Vertex, Connectivity
    
    def gen_vtpxml(self):
        points = self.v
        connectivity = self.p
        Npoints = len(points)
        Npanel = len(connectivity)
        
        print('Converting the mesh file into a Mayavi compatible format xml')
        #create the vtk xml file
        f_n = os.path.join(self.path, '__Mesh.vtp')
        fid = open(f_n, 'w')
        fid.write('<?xml version="1.0"?>\n')
        fid.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        fid.write('  <PolyData>\n')
        fid.write('    <Piece NumberOfPoints="{}" NumberOfVerts="0" NumberOfLines="0"\n'.format(Npoints))
        fid.write('           NumberOfStrips="0" NumberOfPolys="{}">\n'.format(Npanel))
        
        #definition of the vertexes
        fid.write('      <Points>\n')
        fid.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for point in points:
            fid.write('          {} {} {}\n'.format(point[0],point[1],point[2]))
        fid.write('        </DataArray>\n')
        fid.write('      </Points>\n')
        
        #definition of the scalar of each vertex
        fid.write('\n')
        fid.write('      <PointData Scalars="my_scalars">\n')
        fid.write('        <DataArray type="Float32" Name="my_scalars" format="ascii">\n')
        for point in range(Npoints):
            fid.write(' {}'.format(points[point][2]))
        fid.write('\n        </DataArray>\n')
        fid.write('      </PointData>\n')
        
        #definition of the normals
        fid.write('\n')
        fid.write('      <CellData Scalars="cell_scalars" Normals="cell_normals">\n')
        fid.write('        <DataArray type="Int32" Name="cell_scalars" format="ascii">\n')
        for panel in connectivity:
            fid.write('          {} '.format(0))
        fid.write('\n')
        fid.write('        </DataArray>\n')
        
        fid.write('        <DataArray type="Float32" Name="cell_normals"\n')
        fid.write('                   NumberOfComponents="3" format="ascii">\n')
        for el in range(Npanel):
            fid.write('          {} {} {}\n'.format(0,0,1))
        fid.write('        </DataArray>\n')
        fid.write('      </CellData>\n')
        
        
        #definition of the polygons--connectivity and vertex offset
        fid.write('      <Polys>\n')
        fid.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for panel in connectivity:
            fid.write('          {} {} {} {}\n'.format(int(panel[0]),int(panel[1]),int(panel[2]),int(panel[3])))
        fid.write('        </DataArray>\n')
        fid.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        off = 0
        for offset in range(Npanel):
            off+=4
            fid.write('          {}'.format(off))
        fid.write('\n        </DataArray>\n')
        fid.write('      </Polys>\n')
        fid.write('    </Piece>\n')
        fid.write('  </PolyData>\n')
        fid.write('</VTKFile>')
        
        fid.close()
        
        print("vtk xml file (__Mesh.vtp) generated.")
        return 1
        
if __name__ == "__main__":
    fn = "mesh_test.gdf"
    path = ""
    ms = Mesh(fn, path)
    ms.generate_visualiser()
    
