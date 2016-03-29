#version 3.6;

global_settings{ assumed_gamma 1.0 }
global_settings { ambient_light rgb<0, 0, 0> }

#include "colors.inc"
#include "textures.inc"
#include "golds.inc"
#include "lab2_data.inc"

// camera ------------------------------------------------------------------
#declare Cam =camera {ultra_wide_angle angle 45
                                       location  <0.0 , -4.0 ,-30.0>
                                       right     x*image_width/image_height
                                       look_at   <0.0 , 0.0 , 0.0>}
camera{Cam}
// sun ---------------------------------------------------------------------
light_source{ <0,0,0> color rgb <1,1,1>   
              spotlight
              point_at<-1,-1,1>
              translate< 10, 10, -10>
            }
//==========================================================================
#for (vertexIndex, 0, truncated_icosahedron_vertices_count-1)
sphere {
    truncated_icosahedron_vertices[vertexIndex],
    truncated_icosahedron_sphere_radius
    texture {T_Gold_2A}
}
#end
#for (edgeIndex, 0, truncated_icosahedron_edges_count-1)
cylinder {
    truncated_icosahedron_vertices[truncated_icosahedron_edges[edgeIndex].x],
    truncated_icosahedron_vertices[truncated_icosahedron_edges[edgeIndex].y],
    truncated_icosahedron_cylinder_radius
    texture {T_Gold_2A}
}
#end
//==========================================================================
#for (vertexIndex, 0, torus_vertices_count-1)
sphere {
    torus_vertices[vertexIndex],
    torus_sphere_radius
    texture {Jade}
}
#end
#for (edgeIndex, 0, torus_edges_count-1)
cylinder {
    torus_vertices[torus_edges[edgeIndex].x],
    torus_vertices[torus_edges[edgeIndex].y],
    torus_cylinder_radius
    texture {Jade}
}
#end
//--------------------------------- end ------------------------------------
