// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>



#include "src/Triangle.h"
#include "src/Mesh.h"
#include "src/Camera.h"
#include "src/Simplifier.h"
#include "src/Compressor.h"

int weight_type;

enum DisplayMode{ WIRE=0, SOLID=1, LIGHTED_WIRE=2, LIGHTED=3 };

//Transformation made of a rotation and translation
struct Transformation {
    Mat3 rotation;
    Vec3 translation;
};

//Basis ( origin, i, j ,k )
struct Basis {
    inline Basis ( Vec3 const & i_origin,  Vec3 const & i_i, Vec3 const & i_j, Vec3 const & i_k) {
        origin = i_origin; i = i_i ; j = i_j ; k = i_k;
    }

    inline Basis ( ) {
        origin = Vec3(0., 0., 0.);
        i = Vec3(1., 0., 0.) ; j = Vec3(0., 1., 0.) ; k = Vec3(0., 0., 1.);
    }
    Vec3 operator [] (unsigned int ib) {
        if(ib==0) return i;
        if(ib==1) return j;
        return k;}

    Vec3 origin;
    Vec3 i;
    Vec3 j;
    Vec3 k;
};

//Fonction à completer
void collect_one_ring (std::vector<Vec3> const & i_vertices,
                       std::vector< Triangle > const & i_triangles,
                       std::vector<std::vector<unsigned int> > & o_one_ring) {//one-ring of each vertex, i.e. a list of vertices with which it shares an edge
    //Initialiser le vecetur de o_one_ring de la taille du vecteur vertices
    //Parcourir les triangles et ajouter les voisins dans le 1-voisinage
    //Attention verifier que l'indice n'est pas deja present
    //Tous les points opposés dans le triangle sont reliés

    // for(size_t i = 0; i<i_vertices.size(); i++)
    // {
    // 	o_one_ring.push_back(std::vector<unsigned int>());
    // }

    o_one_ring.resize(i_vertices.size());

    for(Triangle triangle : i_triangles)
    {
    	
    	for(int i =0; i<3; i++)
    	{
    		if(std::find(o_one_ring[triangle[i]].begin(), o_one_ring[triangle[i]].end(), i==0? triangle[i+1] : triangle[i-1])== o_one_ring[triangle[i]].end()) 
    		{
    			o_one_ring[triangle[i]].push_back(i==0? triangle[i+1] : triangle[i-1]);
    		}
    		if(std::find(o_one_ring[triangle[i]].begin(), o_one_ring[triangle[i]].end(), i==1? triangle[i+1] : triangle[abs(i-2)]) == o_one_ring[triangle[i]].end())
    		{
    			o_one_ring[triangle[i]].push_back(i==1? triangle[i+1] : triangle[abs(i-2)]);
    		}
    	}

    }

}

//Fonction à completer
void compute_vertex_valences (const std::vector<Vec3> & i_vertices,
                              const std::vector< Triangle > & i_triangles,
                              std::vector<unsigned int> & o_valences ) {
    //Utiliser la fonction collect_one_ring pour récuperer le 1-voisinage
	std::vector<std::vector<unsigned int> > o_one_ring;
	collect_one_ring(i_vertices, i_triangles, o_one_ring);
	for(size_t i = 0; i < i_vertices.size(); i++)
	{
		o_valences.push_back(o_one_ring.at(i).size());
	}

}

//Input mesh loaded at the launch of the application
Mesh mesh;
std::vector< float > mesh_valence_field; //normalized valence of each vertex
Mesh result; //result of simplification;
bool simplified = false;
bool quant = false;

Basis basis;

bool display_normals;
bool display_smooth_normals;
bool display_mesh;
bool display_basis;
bool display_boardingBox;
DisplayMode displayMode;


// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 1600;
static unsigned int SCREENHEIGHT = 900;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;

// ------------------------------------
// File I/O
// ------------------------------------
bool saveOFF( const std::string & filename ,
              std::vector< Vec3 > const & i_vertices ,
              std::vector< Vec3 > const & i_normals ,
              std::vector< Triangle > const & i_triangles,
              std::vector< Vec3 > const & i_triangle_normals ,
              bool save_normals = true ) {
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open()) {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl ;

    unsigned int n_vertices = i_vertices.size() , n_triangles = i_triangles.size();
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v ) {
        myfile << i_vertices[v][0] << " " << i_vertices[v][1] << " " << i_vertices[v][2] << " ";
        if (save_normals) myfile << i_normals[v][0] << " " << i_normals[v][1] << " " << i_normals[v][2] << std::endl;
        else myfile << std::endl;
    }
    for( unsigned int f = 0 ; f < n_triangles ; ++f ) {
        myfile << 3 << " " << i_triangles[f][0] << " " << i_triangles[f][1] << " " << i_triangles[f][2]<< " ";
        if (save_normals) myfile << i_triangle_normals[f][0] << " " << i_triangle_normals[f][1] << " " << i_triangle_normals[f][2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}

void openOFF( std::string const & filename,
              std::vector<Vec3> & o_vertices,
              std::vector<Vec3> & o_normals,
              std::vector< Triangle > & o_triangles,
              std::vector< Vec3 > & o_triangle_normals,
              bool load_normals = true )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        exit(1);
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    o_vertices.clear();
    o_normals.clear();

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        float x , y , z ;

        myfile >> x >> y >> z ;
        o_vertices.push_back( Vec3( x , y , z ) );

        if( load_normals ) {
            myfile >> x >> y >> z;
            o_normals.push_back( Vec3( x , y , z ) );
        }
    }

    o_triangles.clear();
    o_triangle_normals.clear();
    for( int f = 0 ; f < n_faces ; ++f )
    {
        int n_vertices_on_face;
        myfile >> n_vertices_on_face;

        if( n_vertices_on_face == 3 )
        {
            unsigned int _v1 , _v2 , _v3;
            myfile >> _v1 >> _v2 >> _v3;

            o_triangles.push_back(Triangle( _v1, _v2, _v3 ));

            if( load_normals ) {
                float x , y , z ;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back( Vec3( x , y , z ) );
            }
        }
        else if( n_vertices_on_face == 4 )
        {
            unsigned int _v1 , _v2 , _v3 , _v4;
            myfile >> _v1 >> _v2 >> _v3 >> _v4;

            o_triangles.push_back(Triangle(_v1, _v2, _v3 ));
            o_triangles.push_back(Triangle(_v1, _v3, _v4));
            if( load_normals ) {
                float x , y , z ;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back( Vec3( x , y , z ) );
            }

        }
        else
        {
            std::cout << "We handle ONLY *.off files with 3 or 4 vertices per face" << std::endl;
            myfile.close();
            exit(1);
        }
    }

}

// ------------------------------------
// Cho's Method
// ------------------------------------

void normalize(std::vector<Vec3> &i_values, float min, float max)
{
    for (int i = 0; i < i_values.size(); ++i)
    {
        i_values[i][2] = (i_values[i][2]-min)/(max - min);
    }
}

void unnormalize(std::vector<Vec3> &normalizedValues, float min, float max)
{
    for (int i = 0; i < normalizedValues.size(); ++i)
    {
        normalizedValues[i][2] = normalizedValues[i][2] * (max - min) + min;
    }
}

std::vector<std::vector<int>> computeBins(std::vector<Vec3> &values, float min, float max, int N)
{
    std::vector<std::vector<int>> bins = std::vector<std::vector<int>>(N);
    int L = values.size();
    for (int i = 0; i < L; ++i)
    {
        float value = values[i][2];
        int binIndex = -1;
        for (int n = 0; n < N; ++n)
        {
            float binMin = min + ((max - min)/N * n);
            float binMax = min + ((max - min)/N * (n+1));
            if(value >= binMin && value < binMax)
            {
                binIndex = n;
                break;
            }
        }
        if(binIndex >= 0) bins[binIndex].push_back(i);
    }
    return bins;
}

void embed(int omega, std::vector<std::vector<int>> &bins, int n, float alpha)
{
    float kn = 1;
    std::vector<Vec3> transformedBin(bins[n].size());
    float rhonj_;
    for (int j = 0; j < bins[n].size(); ++j)
    {
        transformedBin[j] = mesh.vertices[bins[n][j]];
        rhonj_ = pow(mesh.vertices[bins[n][j]][2] , kn);
        transformedBin[j][2] = rhonj_;
    }

    float mean = 0.0;

    for (int j = 0; j < transformedBin.size(); ++j)
    {
        mean += transformedBin[j][2];
    }
    mean /= (float)transformedBin.size();

    while ((omega == 1) ? (mean < (0.5f + alpha)) : (mean > (0.5f - alpha)))
    {
        kn += (omega == 1) ? -0.01f : 0.01f;

        for (int j = 0; j < transformedBin.size(); ++j)
        {
            rhonj_ = pow(transformedBin[j][2], kn);
            transformedBin[j][2]= rhonj_;
        }

        mean = 0.0;

        for (int j = 0; j < transformedBin.size(); ++j)
        {
            mean += transformedBin[j][2];
        }
        mean /= (float)transformedBin.size();

        
    }

    for (int i = 0; i < bins[n].size(); ++i)
    {
        mesh.vertices[bins[n][i]][2] = transformedBin[i][2];
    }

}


void watermark(Mesh &mesh, std::vector<int> &message, int nb_bins, float alpha)
{
    // Convert to spherical and store norms
    std::vector<float> norms(mesh.vertices.size());
    std::vector<Vec3> sphericals(mesh.vertices.size());
    std::vector<float> minmax = mesh.radialMinMax();

    for (int i = 0; i < mesh.vertices.size(); ++i)
    {
        mesh.vertices[i] = Vec3::EuclideanCoordinatesToSpherical(mesh.vertices[i]);
    }
    std::vector<std::vector<int>> bins = computeBins(mesh.vertices, minmax[0], minmax[1], nb_bins);
    norms.clear();
    for (int n = 0; n < bins.size(); ++n)
    {
        // compute min and max of each bin
        printf("Bin %d\n", n);
        printf("Computing min max\n");
        float min, max;
        min = FLT_MAX;
        max = -FLT_MAX;
        for (int j = 0; j < bins[n].size(); ++j)
        {
            if(mesh.vertices[bins[n][j]][2] < min) min = mesh.vertices[bins[n][j]][2] ;
            if(mesh.vertices[bins[n][j]][2]  > max) max = mesh.vertices[bins[n][j]][2] ;
        }

        //embeding
        printf("embeding\n");
        normalize(mesh.vertices, min, max);
        embed(message[n], bins, n, alpha);
        unnormalize(mesh.vertices, min, max);
        
        
    }
    for (int i = 0; i < mesh.vertices.size(); ++i)
    {
        /* code */
        mesh.vertices[i] = Vec3::SphericalCoordinatesToEuclidean(mesh.vertices[i]);
    } 

}

void extraction(Mesh &mesh, int nb_bins, float alpha)
{
    std::vector<int> message = std::vector<int>();
    std::vector<float> minmax = mesh.radialMinMax();

    for (int i = 0; i < mesh.vertices.size(); ++i)
    {
        mesh.vertices[i] = Vec3::EuclideanCoordinatesToSpherical(mesh.vertices[i]);
    }
    std::vector<std::vector<int>> bins = computeBins(mesh.vertices, minmax[0], minmax[1], nb_bins);
    for (int n = 0; n < bins.size(); ++n)
    {
        // compute min and max of each bin
        printf("Bin %d\n", n);
        printf("Computing min max\n");
        float min, max;
        min = FLT_MAX;
        max = -FLT_MAX;
        for (int j = 0; j < bins[n].size(); ++j)
        {
            if(mesh.vertices[bins[n][j]][2] < min) min = mesh.vertices[bins[n][j]][2] ;
            if(mesh.vertices[bins[n][j]][2]  > max) max = mesh.vertices[bins[n][j]][2] ;
        }

        
        normalize(mesh.vertices, min, max);
        float mean = 0.0f;
        for (int i = 0; i < bins[n].size(); ++i)
        {
            mean += mesh.vertices[bins[n][i]][2];
        }
        mean /= bins[n].size();
        if(mean > 0.5f) message.push_back(1);
        else message.push_back(-1);
        unnormalize(mesh.vertices, min, max);
        for (int i = 0; i < mesh.vertices.size(); ++i)
        {
            /* code */
            mesh.vertices[i] = Vec3::SphericalCoordinatesToEuclidean(mesh.vertices[i]);
        } 


            
        
    }
    for (int i = 0; i < message.size(); ++i)
    {
        printf("%d, ", message[i]);
    }
    printf("\n");
}

// void testNormalize()
// {
//     std::vector<float> values = {2.5, 3.4, 5.6, 1.2};
//     values = normalize(values, 1.2, 5.6);
//     for (int i = 0; i < values.size(); ++i)
//     {
//         printf("%f,", values[i]);
//     }
//     printf("\nnormalize ok \n");
// }
// void testUnormalize()
// {
//     std::vector<float> values = {0.295455,0.500000,1.000000,0.000000};
//     values = unnormalize(values, 1.2, 5.6);
//     for (int i = 0; i < values.size(); ++i)
//     {
//         printf("%f,", values[i]);
//     }
//     printf("\nunnormalize ok \n");
// }
// void testComputeBins()
// {
//     std::vector<float> values = {2.5, 3.4, 5.6, 1.2, 2.2, 3.3, 4.4, 5.4, 2.345, 1.234, 4.567};
//     std::vector<std::vector<float>> bins = computeBins(values, 1.2, 5.6, 3);
//     for (int i = 0; i < bins.size(); ++i)
//     {
//         printf("bin %d\n", i);
//         for (int j = 0; j < bins[i].size(); ++j)
//         {
//             printf("---  %f\n", bins[i][j]);
//         }
//     }
    

// }
// void testEmbed()
// {
//     std::vector<float> values = {2.5, 3.4, 5.6, 1.2, 2.2, 3.3, 4.4, 5.4, 2.345, 1.234, 4.567};
//     std::vector<std::vector<float>> bins = computeBins(values, 1.2, 5.6, 3);
//     float min, max;
//     min = FLT_MAX;
//     max = -FLT_MAX;
//     for (int j = 0; j < bins[0].size(); ++j)
//     {
//         if(bins[0][j] < min) min = bins[0][j];
//         if(bins[0][j] > max) max = bins[0][j];
//     }
//     printf("min %f, max %f\n", min, max);
//     bins[0] = normalize(bins[0], min, max);
//     embed(-1, bins, 0, 0.1);
//     bins[0] = unnormalize(bins[0], min, max);
//     for (int i = 0; i < bins[0].size(); ++i)
//     {
//         printf("%f\n", bins[0][i]);
//     }

// }

// ------------------------------------
// Application initialization
// ------------------------------------
void initLight () {
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    glCullFace (GL_BACK);
    glDisable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    display_normals = false;
    display_mesh = true;
    display_smooth_normals = true;
    displayMode = LIGHTED;
    display_basis = false;
    display_boardingBox = false;
}


// ------------------------------------
// Rendering.
// ------------------------------------

void drawVector( Vec3 const & i_from, Vec3 const & i_to ) {

    glBegin(GL_LINES);
    glVertex3f( i_from[0] , i_from[1] , i_from[2] );
    glVertex3f( i_to[0] , i_to[1] , i_to[2] );
    glEnd();
}

void drawAxis( Vec3 const & i_origin, Vec3 const & i_direction ) {

    glLineWidth(4); // for example...
    drawVector(i_origin, i_origin + i_direction);
}

void drawReferenceFrame( Vec3 const & origin, Vec3 const & i, Vec3 const & j, Vec3 const & k ) {

    glDisable(GL_LIGHTING);
    glColor3f( 0.8, 0.2, 0.2 );
    drawAxis( origin, i );
    glColor3f( 0.2, 0.8, 0.2 );
    drawAxis( origin, j );
    glColor3f( 0.2, 0.2, 0.8 );
    drawAxis( origin, k );
    glEnable(GL_LIGHTING);

}

void drawReferenceFrame( Basis & i_basis ) {
    drawReferenceFrame( i_basis.origin, i_basis.i, i_basis.j, i_basis.k );
}

typedef struct {
    float r;       // ∈ [0, 1]
    float g;       // ∈ [0, 1]
    float b;       // ∈ [0, 1]
} RGB;



RGB scalarToRGB( float scalar_value ) //Scalar_value ∈ [0, 1]
{
    RGB rgb;
    float H = scalar_value*360., S = 1., V = 0.85,
            P, Q, T,
            fract;

    (H == 360.)?(H = 0.):(H /= 60.);
    fract = H - floor(H);

    P = V*(1. - S);
    Q = V*(1. - S*fract);
    T = V*(1. - S*(1. - fract));

    if      (0. <= H && H < 1.)
        rgb = (RGB){.r = V, .g = T, .b = P};
    else if (1. <= H && H < 2.)
        rgb = (RGB){.r = Q, .g = V, .b = P};
    else if (2. <= H && H < 3.)
        rgb = (RGB){.r = P, .g = V, .b = T};
    else if (3. <= H && H < 4.)
        rgb = (RGB){.r = P, .g = Q, .b = V};
    else if (4. <= H && H < 5.)
        rgb = (RGB){.r = T, .g = P, .b = V};
    else if (5. <= H && H < 6.)
        rgb = (RGB){.r = V, .g = P, .b = Q};
    else
        rgb = (RGB){.r = 0., .g = 0., .b = 0.};

    return rgb;
}

void drawSmoothTriangleMesh( Mesh const & i_mesh , bool draw_field = false ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_mesh.triangles.size(); ++tIt) {

        for(unsigned int i = 0 ; i < 3 ; i++) {
            const Vec3 & p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; //Vertex position
            const Vec3 & n = i_mesh.normals[i_mesh.triangles[tIt][i]]; //Vertex normal

            if( draw_field && mesh_valence_field.size() > 0 ){
                RGB color = scalarToRGB( mesh_valence_field[i_mesh.triangles[tIt][i]] );
                glColor3f( color.r, color.g, color.b );
            }
            glNormal3f( n[0] , n[1] , n[2] );
            glVertex3f( p[0] , p[1] , p[2] );
        }
    }
    glEnd();

}

void drawTriangleMesh( Mesh const & i_mesh , bool draw_field = false  ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_mesh.triangles.size(); ++tIt) {
        const Vec3 & n = i_mesh.triangle_normals[ tIt ]; //Triangle normal
        for(unsigned int i = 0 ; i < 3 ; i++) {
            const Vec3 & p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; //Vertex position

            if( draw_field ){
                RGB color = scalarToRGB( mesh_valence_field[i_mesh.triangles[tIt][i]] );
                glColor3f( color.r, color.g, color.b );
            }
            glNormal3f( n[0] , n[1] , n[2] );
            glVertex3f( p[0] , p[1] , p[2] );
        }
    }
    glEnd();

}

void drawMesh( Mesh const & i_mesh , bool draw_field = false ){
    if(display_smooth_normals)
        drawSmoothTriangleMesh(i_mesh, draw_field) ; //Smooth display with vertices normals
    else
        drawTriangleMesh(i_mesh, draw_field) ; //Display with face normals
}

void drawVectorField( std::vector<Vec3> const & i_positions, std::vector<Vec3> const & i_directions ) {
    glLineWidth(1.);
    for(unsigned int pIt = 0 ; pIt < i_directions.size() ; ++pIt) {
        Vec3 to = i_positions[pIt] + 0.02*i_directions[pIt];
        drawVector(i_positions[pIt], to);
    }
}

void drawNormals(Mesh const& i_mesh){

    if(display_smooth_normals){
        drawVectorField( i_mesh.vertices, i_mesh.normals );
    } else {
        std::vector<Vec3> triangle_baricenters;
        for ( const Triangle& triangle : i_mesh.triangles ){
            Vec3 triangle_baricenter (0.,0.,0.);
            for( unsigned int i = 0 ; i < 3 ; i++ )
                triangle_baricenter += i_mesh.vertices[triangle[i]];
            triangle_baricenter /= 3.;
            triangle_baricenters.push_back(triangle_baricenter);
        }

        drawVectorField( triangle_baricenters, i_mesh.triangle_normals );
    }
}

void drawBoardingBox(Mesh &i_mesh)
{
	glPointSize(20.);
	glBegin(GL_POINTS);
	std::vector<Vec3> boardingbox = i_mesh.boardingbox;
	for(size_t i = 0; i<boardingbox.size(); i++)
	{
		glVertex3f( boardingbox[i][0] , boardingbox[i][1] , boardingbox[i][2] );
	}
	glEnd();

}

//Draw fonction
void draw () {



    if(displayMode == LIGHTED || displayMode == LIGHTED_WIRE){

        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_LIGHTING);

    }  else if(displayMode == WIRE){

        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glDisable (GL_LIGHTING);

    }  else if(displayMode == SOLID ){
        glDisable (GL_LIGHTING);
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    }

    glColor3f(0.8,1,0.8);
    if(simplified) drawMesh(result, true);
    else drawMesh(mesh, true);

    if(displayMode == SOLID || displayMode == LIGHTED_WIRE){
        glEnable (GL_POLYGON_OFFSET_LINE);
        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glLineWidth (1.0f);
        glPolygonOffset (-2.0, 1.0);

        glColor3f(0.,0.,0.);
        drawMesh(mesh, false);

        glDisable (GL_POLYGON_OFFSET_LINE);
        glEnable (GL_LIGHTING);
    }

    if(display_boardingBox)
    {
    	glColor3f(0.5, 0.5, 0.);
    	drawBoardingBox(mesh);
    }

    if(quant) 
    {
        std::vector<Vec3> init_v = mesh.vertices;
        // for (int qp = 5; qp <= 30; ++qp)
        // {
        //     mesh.vertices = Compressor::quant_dequant(mesh, qp);
        //     double rmse = Compressor::RMSE(init_v, mesh.vertices);
        //     printf("%d %f\n", qp, rmse);
        //     mesh.vertices = init_v;
        // }
        mesh.vertices = Compressor::quant_dequant(mesh, 10);
        quant = false;
    }
    //Compressor::dequantify(mesh, 30);


    glDisable(GL_LIGHTING);
    if(display_normals){
        glColor3f(1.,0.,0.);
        drawNormals(mesh);
    }

    if( display_basis ){
        drawReferenceFrame(basis);
    }
    glEnable(GL_LIGHTING);


}

void changeDisplayMode(){
    if(displayMode == LIGHTED)
        displayMode = LIGHTED_WIRE;
    else if(displayMode == LIGHTED_WIRE)
        displayMode = SOLID;
    else if(displayMode == SOLID)
        displayMode = WIRE;
    else
        displayMode = LIGHTED;
}

void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    glutPostRedisplay ();
}

// ------------------------------------
// User inputs
// ------------------------------------
//Keyboard event
void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;


    case 'w': //Change le mode d'affichage
        changeDisplayMode();
        break;


    case 'b': //Toggle basis display
        display_basis = !display_basis;
        break;

    case 'n': //Press n key to display normals
        display_normals = !display_normals;
        break;

    case '1': //Toggle loaded mesh display
        display_mesh = !display_mesh;
        break;

    case 'v': //Switches between face normals and vertices normals
        display_smooth_normals = !display_smooth_normals;
        break;

    case 's':
        Simplifier::simplify(64,mesh, result);
        simplified = true;
        break;

    case 'x':
    	//display_boardingBox = !display_boardingBox;
        quant = !quant;
    	break;

    case '+': //Changes weight type: 0 uniforme, 1 aire des triangles, 2 angle du triangle
        weight_type ++;
        if(weight_type == 3) weight_type = 0;
        mesh.computeVerticesNormals(weight_type); //recalcul des normales avec le type de poids choisi
        break;

    default:
        break;
    }
    idle ();
}

//Mouse events
void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }

    idle ();
}

//Mouse motion, update camera
void motion (int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize (w, h);
}



// ------------------------------------
// Start of graphical application
// ------------------------------------
int main (int argc, char ** argv) {
    if (argc > 2) {
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("TP HAI918I");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);

    //Mesh loaded with precomputed normals
    openOFF("data/bunny.off", mesh.vertices, mesh.normals, mesh.triangles, mesh.triangle_normals, false);

    //Completer les fonction de calcul de normals
    mesh.computeNormals(weight_type);
    std::vector<int> message = {1, -1, 1, -1, -1, -1, -1, 1, 1};
    watermark(mesh, message, 9, 0.49f);
    extraction(mesh, 9, 0.49f);
    

    // Tests
    // testNormalize();
    // testUnormalize();
    // testComputeBins();
    // testEmbed();
    

    glutMainLoop ();
    return EXIT_SUCCESS;
}

