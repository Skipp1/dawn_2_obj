#include "dawn_2_obj.h"
#include <memory>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

std::array<double, 3> array_add(const std::array<double, 3> &a, const std::array<double, 3> b) {
	return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

size_t model::add_unique_vertex(std::vector<std::array<double, 3>> &v, const std::array<double, 3> &p) {
	auto it = std::find(v.begin(), v.end(), p);
	if (it == v.end()) {
		v.push_back(p);
		return v.size();
	} else {
		return std::distance(v.begin(), it)+1;
	}
}


void model::add_box(const std::array<double, 3> &s) {
	std::stringstream ss;
	ss << "g group_" << group_no++ << std::endl;
	ss << "usemtl colour_" << add_unique_vertex(colour_db, rgb) << std::endl;
	ss << "f " << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]-s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]-s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]+s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]+s[1], o[2]-s[2]})
	   << std::endl;
	ss << "f " << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]-s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]-s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]+s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]+s[1], o[2]+s[2]})
	   << std::endl;
	ss << "f " << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]-s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]+s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]+s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]-s[1], o[2]+s[2]})
	   << std::endl;
	ss << "f " << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]-s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]+s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]+s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]-s[1], o[2]+s[2]})
	   << std::endl;
	ss << "f " << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]-s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]-s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]-s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]-s[1], o[2]-s[2]})
	   << std::endl;
	ss << "f " << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]+s[1], o[2]-s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]-s[0], o[1]+s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]+s[1], o[2]+s[2]})
	   << " "  << add_unique_vertex(vertex_db, {o[0]+s[0], o[1]+s[1], o[2]-s[2]})
	   << std::endl;
	obj_db.push_back(ss.str());
}


void model::add_line(const std::vector<std::array<double, 3>> &line) {
	std::stringstream ss;
	ss << "g group_" << group_no++ << std::endl;
	ss << "usemtl colour_" << add_unique_vertex(colour_db, rgb) << std::endl;
	ss << "l";
	for (const auto &p : line ) {
		ss << " " << add_unique_vertex(vertex_db, array_add(p, o));
	}
	ss << std::endl;
	obj_db.push_back(ss.str());
}


void model::write(const std::string &filename) {
	std::ofstream fp_obj(filename + ".obj");
	fp_obj << std::setprecision(8) << std::fixed;
	for (const auto &v : vertex_db) {
		fp_obj << "v " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	for (const auto &s : obj_db ) { 
		fp_obj << s << std::endl;
	}
	fp_obj.close();
	
	std::ofstream fp_mat(filename + ".mtl");
	fp_mat << std::setprecision(8) << std::fixed;
	for ( size_t i=1; const auto &v : colour_db ) { 
		fp_mat << "newmtl color_" << i++ << std::endl;
		fp_mat << "Kd " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	fp_mat.close();
}


extern "C" void *construct_model() {
	return static_cast<void *>(new model);
}

extern "C" void set_rgb(void * _self, double r, double g, double b){
	auto self = static_cast<model *>(_self); 
	self->rgb = {r, g, b};
}

extern "C" void set_origin(void * _self, double x, double y, double z){
	auto self = static_cast<model *>(_self); 
	self->o = {x, y, z};
}

extern "C" void set_vertex(void * _line, double x, double y, double z){
	auto line = static_cast<std::vector<std::array<double, 3>> *>(_line); 
	line->push_back({x, y, z});
}

extern "C" void *polyline() {
	return static_cast<void *>(new std::vector<std::array<double, 3>>);
}

extern "C" void add_line(void *_self, void *_line) {
	auto self = static_cast<model *>(_self); 
	auto line = static_cast<std::vector<std::array<double, 3>> *>(_line); 
	self->add_line(*line);
	delete line;
}

extern "C" void add_box(void *_self, double dx, double dy, double dz) {
	auto self = static_cast<model *>(_self); 
	self->add_box({dz, dy, dz});
}

extern "C" void write_obj(void *_self, const char *filename) {
	auto self = static_cast<model *>(_self); 
	self->write(filename);
}
