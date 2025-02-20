#ifdef __cplusplus
#include <string> 
#include <array> 
#include <vector>

struct polyhedron_t {
	std::vector<std::array<double, 3>> v;
	std::vector<std::array<long int, 3>> f;
};

class model {
	public:
	model(std::string fn):filename(fn){}
		
	void add_line(const std::vector<std::array<double, 3>> &line);
	void add_box(const std::array<double, 3> &oset, const std::array<double, 3> &s);
	void add_sphere( const std::array<double, 3> &o, double r, double res);
	void add_polyhedron(const polyhedron_t &p);
	void add_basis(const std::array<double, 3> &e1, const std::array<double, 3> &e2);
	void write(const std::string &filename);
	
	std::array<double, 3> rgb;
	std::array<double, 3> origin;
	
	const std::string filename;
	double marker_res = 10;
	double model_res = 10;
	
	private:
	size_t add_unique_vertex(std::vector<std::array<double, 3>> &v, const std::array<double, 3> &p);
	
	std::vector<std::array<double, 3>> vertex_db, colour_db; 
	std::vector<std::string> obj_db;
	
	std::array<std::array<double, 3>, 3> xform = { std::array<double, 3>{1, 0, 0}
	                                             , std::array<double, 3>{0, 1, 0}
	                                             , std::array<double, 3>{0, 0, 1}};
	
	size_t group_no = 0;
};

extern "C" {
#endif 

void yytext_dup( const char *c );
	
void *construct_model(int argc, char **argv);

void set_rgb(void *self, double r, double g, double b);
void set_origin(void *self, double x, double y, double z);
void set_basis(void *_self, double e1x, double e1y, double e1z, double e2x, double e2y, double e2z);
	
void *polyline();
void line_add_vertex(void *line, double x, double y, double z); 
void add_line(void *self, void *line );

void add_box(void *self, double dx, double dy, double dz);
void add_box_mark(void *self, double x, double y, double z, double r);

void add_sphere(void *self, double x, double y, double z, double r);
void add_sphere_mark(void *self, double x, double y, double z, double r);

void *polyhedron();
void polyhedron_add_vertex(void *polyhedron, double x, double y, double z);
void polyhedron_add_face(void *polyhedron, long int i, long int j, long int k);
void add_polyhedron(void *self, void *polyhedron);

void write_obj(void *self);
#ifdef __cplusplus
}
#endif
