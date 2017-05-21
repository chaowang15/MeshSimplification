/*
Fast and Robust Detection of Crest Lineson Meshes C++ code
Copyright:(c) Shin Yoshizawa, 2004
E-mail: shin.yoshizawa@mpi-sb.mpg.de
URL: http://www.mpi-sb.mpg.de/~shin
Affiliation: Max-Planck-Institut fuer Informatik: Computer Graphics Group 
 Stuhlsatzenhausweg 85, 66123 Saarbruecken, Germany
 Phone +49 681 9325-408 Fax +49 681 9325-499 

 All right is reserved by Shin Yoshizawa.
This C++ sources are allowed for only primary user of 
research and educational purposes. Don't use secondary: copy, distribution, 
diversion, business purpose, and etc.. 
 */
#include<stdio.h>
#include<math.h>
#include"Point3d.h"
#include"PointTool.h"
#include"IDList.h"
#include"PolarList.h"
#include"IDSet.h"
#include"Polyhedron.h"

int main(int argc,char *argv[]){
  
  Polyhedron *mymesh = new Polyhedron();  
  mymesh->readmesh(argv[1]);
  delete mymesh;
  return 0;

}

