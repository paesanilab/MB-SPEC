#include <cassert>

#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>

#include "read-lammpstrj.h"
#include "io-tools.h"

////////////////////////////////////////////////////////////////////////////////

namespace opt {
    

//----------------------------------------------------------------------------//

void read_lammpstrj_frame(size_t& lineno, size_t& imcon, std::istream& is,
                          opt::lammpstrj_frame& frame, double t,
                          bool require_dipoles)
{
    assert(is && !is.eof());
    assert(frame.natm == 0); // empty

    // Reading a line and throwing it
    std::string line;
    std::getline(is, line); 
    if (is.eof()) // check eof before defining/allocating variables
        return;
    ++lineno;

    // Reading a line and saving step number
    std::getline(is, line); 
    if (is.eof()) 
        return;
    ++lineno;
    std::istringstream iss(line);
    iss >> frame.timestep;
    frame.time = t;

    // Reading a line and throwing it
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

    // Reading a line and saving total atom number
    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;
    iss.clear();
    iss.str(line);
    iss >> frame.natm;

    frame.xyz        = new double[3*frame.natm];
    frame.name       = new char[frame.natm*MAX_CHAR_LENGTH];
    frame.atomic_num = new int[frame.natm];
    frame.mass       = new double[frame.natm];
    frame.cell       = new double[9];
    frame.rcell       = new double[9];
    frame.dipind       = new double[4*frame.natm];
    frame.dipmol       = new double[4*frame.natm];
    frame.molecid       = new int[frame.natm];

    std::fill(frame.name, frame.name + frame.natm*MAX_CHAR_LENGTH, '0');
    std::fill(frame.cell, frame.cell + 9, 0);
    std::fill(frame.rcell, frame.rcell + 9, 0);
    std::fill(frame.dipind, frame.dipind + 4*frame.natm, 0.0);
    std::fill(frame.dipmol, frame.dipmol + 4*frame.natm, 0.0);
    std::fill(frame.molecid, frame.molecid + frame.natm, 0);
    std::fill(frame.atomic_num, frame.atomic_num + frame.natm, 0);

    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;
    if (line.rfind("ITEM: BOX BOUNDS", 0) != 0)
        throw std::runtime_error("expected 'ITEM: BOX BOUNDS'");
    const bool general_triclinic = line.find("abc") != std::string::npos;
    const bool restricted_triclinic = !general_triclinic &&
        line.find("xy") != std::string::npos &&
        line.find("xz") != std::string::npos &&
        line.find("yz") != std::string::npos;

    std::array<std::vector<double>, 3> box_rows;
    for (size_t i = 0; i < 3; ++i) {
        std::getline(is, line);
        if (is.eof())
            return;
        ++lineno;
        iss.clear();
        iss.str(line);
        for (double value = 0.0; iss >> value; )
            box_rows[i].push_back(value);
    }

    std::array<double, 3> origin = {{0.0, 0.0, 0.0}};
    if (imcon < 3) {
        if (restricted_triclinic || general_triclinic ||
            box_rows[0].size() != 2 || box_rows[1].size() != 2 ||
            box_rows[2].size() != 2) {
            throw std::runtime_error("box with imcon = " + std::to_string(imcon) +
                                     " must be orthogonal");
        }
        origin[0] = box_rows[0][0];
        origin[1] = box_rows[1][0];
        origin[2] = box_rows[2][0];
        frame.cell[0] = box_rows[0][1] - box_rows[0][0];
        frame.cell[4] = box_rows[1][1] - box_rows[1][0];
        frame.cell[8] = box_rows[2][1] - box_rows[2][0];
        frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];
        if (frame.vol <= 0.0)
            throw std::runtime_error("orthogonal LAMMPS box must have positive volume");
        frame.rcell[0] = 1.0/frame.cell[0];
        frame.rcell[4] = 1.0/frame.cell[4];
        frame.rcell[8] = 1.0/frame.cell[8];
    } else if (imcon == 3) {
        if (restricted_triclinic) {
            if (box_rows[0].size() != 3 || box_rows[1].size() != 3 ||
                box_rows[2].size() != 3)
                throw std::runtime_error("restricted triclinic LAMMPS box must have three values per row");
            const double xy = box_rows[0][2];
            const double xz = box_rows[1][2];
            const double yz = box_rows[2][2];
            origin[0] = box_rows[0][0] - std::min({0.0, xy, xz, xy + xz});
            origin[1] = box_rows[1][0] - std::min(0.0, yz);
            origin[2] = box_rows[2][0];
            const double xhi = box_rows[0][1] - std::max({0.0, xy, xz, xy + xz});
            const double yhi = box_rows[1][1] - std::max(0.0, yz);
            frame.cell[0] = xhi - origin[0];
            frame.cell[3] = xy;
            frame.cell[4] = yhi - origin[1];
            frame.cell[6] = xz;
            frame.cell[7] = yz;
            frame.cell[8] = box_rows[2][1] - origin[2];
        } else if (general_triclinic) {
            if (box_rows[0].size() != 4 || box_rows[1].size() != 4 ||
                box_rows[2].size() != 4)
                throw std::runtime_error("general triclinic LAMMPS box must have four values per row");
            for (size_t row = 0; row < 3; ++row) {
                for (size_t component = 0; component < 3; ++component)
                    frame.cell[3*row + component] = box_rows[row][component];
                origin[row] = box_rows[row][3];
            }
        } else {
            throw std::runtime_error("imcon = 3 requires a triclinic LAMMPS box");
        }
        const double cross_x = frame.cell[4]*frame.cell[8] - frame.cell[5]*frame.cell[7];
        const double cross_y = -frame.cell[1]*frame.cell[8] + frame.cell[7]*frame.cell[2];
        const double cross_z = frame.cell[1]*frame.cell[5] - frame.cell[4]*frame.cell[2];
        frame.vol = frame.cell[0]*cross_x + frame.cell[3]*cross_y + frame.cell[6]*cross_z;
        if (frame.vol <= 0.0)
            throw std::runtime_error("triclinic LAMMPS box must have positive volume");
        calculate_rcell(lineno, frame.cell, frame.rcell);
    } else {
        throw std::runtime_error("Do not know how to calculate volume for imcon = " +
                                 std::to_string(imcon));
    }

    std::getline(is, line); 
    if (is.eof())
        return;
    ++lineno;

    const std::string atom_prefix = "ITEM: ATOMS";
    if (line.rfind(atom_prefix, 0) != 0)
        throw std::runtime_error("expected 'ITEM: ATOMS'");
    std::unordered_map<std::string, size_t> columns;
    iss.clear();
    iss.str(line.substr(atom_prefix.size()));
    for (std::string column; iss >> column; ) {
        if (columns.count(column) != 0)
            throw std::runtime_error("duplicate LAMMPS atom column '" + column + "'");
        columns[column] = columns.size();
    }

    const auto column_index = [&columns](const std::string& name) {
        const std::unordered_map<std::string, size_t>::const_iterator found =
            columns.find(name);
        return found == columns.end() ? columns.size() : found->second;
    };
    const std::array<std::array<std::string, 3>, 4> coordinate_names = {{{{
        "x", "y", "z"}}, {{"xu", "yu", "zu"}},
        {{"xs", "ys", "zs"}}, {{"xsu", "ysu", "zsu"}}}};
    size_t coordinate_mode = coordinate_names.size();
    for (size_t mode = 0; mode < coordinate_names.size(); ++mode) {
        const bool complete = column_index(coordinate_names[mode][0]) < columns.size() &&
            column_index(coordinate_names[mode][1]) < columns.size() &&
            column_index(coordinate_names[mode][2]) < columns.size();
        if (!complete)
            continue;
        if (coordinate_mode != coordinate_names.size())
            throw std::runtime_error("LAMMPS atom header must provide exactly one complete coordinate triplet");
        coordinate_mode = mode;
    }
    if (coordinate_mode == coordinate_names.size())
        throw std::runtime_error("LAMMPS atom header must provide x/y/z, xu/yu/zu, xs/ys/zs, or xsu/ysu/zsu");

    std::array<size_t, 3> coordinate_columns;
    for (size_t component = 0; component < 3; ++component)
        coordinate_columns[component] = column_index(coordinate_names[coordinate_mode][component]);
    const bool scaled_coordinates = coordinate_mode == 2 || coordinate_mode == 3;
    const std::array<std::string, 6> dipole_names = {{
        "v_dipx_perm", "v_dipy_perm", "v_dipz_perm",
        "v_dipx_ind", "v_dipy_ind", "v_dipz_ind"}};
    std::array<size_t, 6> dipole_columns;
    for (size_t component = 0; component < dipole_names.size(); ++component)
        dipole_columns[component] = column_index(dipole_names[component]);
    if (require_dipoles && std::any_of(dipole_columns.begin(), dipole_columns.end(),
                                       [&columns](size_t column) { return column >= columns.size(); })) {
        throw std::runtime_error("LAMMPS atom header must provide permanent and induced dipole components");
    }
    const size_t molecule_column = column_index("mol");
    const size_t type_column = column_index("type");

    for (size_t n = 0; n < frame.natm; ++n) {
        if (is.eof()) 
    	    opt::print_error(lineno, 1);

        std::getline(is, line);
        ++lineno;

        std::vector<double> values;
        iss.clear();
        iss.str(line);
        for (double value = 0.0; iss >> value; )
            values.push_back(value);
        if (values.size() != columns.size())
            throw std::runtime_error("LAMMPS atom record at line " +
                                     std::to_string(lineno) + " does not match its header");

        double coordinate[3];
        for (size_t component = 0; component < 3; ++component)
            coordinate[component] = values[coordinate_columns[component]];
        if (scaled_coordinates) {
            for (size_t component = 0; component < 3; ++component) {
                frame.xyz[3*n + component] = origin[component] +
                    frame.cell[component]*coordinate[0] +
                    frame.cell[3 + component]*coordinate[1] +
                    frame.cell[6 + component]*coordinate[2];
            }
        } else {
            for (size_t component = 0; component < 3; ++component)
                frame.xyz[3*n + component] = coordinate[component];
        }
        if (molecule_column < columns.size())
            frame.molecid[n] = static_cast<int>(values[molecule_column]);
        if (type_column < columns.size())
            frame.atomic_num[n] = static_cast<int>(values[type_column]);
        if (require_dipoles) {
            for (size_t component = 0; component < 3; ++component) {
                frame.dipmol[4*n + component] = values[dipole_columns[component]];
                frame.dipind[4*n + component] = values[dipole_columns[3 + component]];
            }
        }
    }
    
    iss.clear();
}

//----------------------------------------------------------------------------//

size_t read_lammpstrj 
	(const char * in, std::vector<opt::lammpstrj_frame>& frames, size_t imcon,
         double dump, double trajectory_duration_ps) {


	std::ifstream ifs(in);
	if (!ifs) {
		std::ostringstream oss;
		oss << "could not open " << in << " for reading\n";
		throw std::runtime_error(oss.str());
	}

	size_t lineno(0);
	size_t frameno(0);

        
	while((!ifs.eof())){
		if (trajectory_duration_ps > 0.0 &&
                    dump*frameno > trajectory_duration_ps + 1.0e-12) break;
		opt::lammpstrj_frame this_frame;
            
		read_lammpstrj_frame(lineno, imcon, ifs, this_frame, dump * frameno);
            
		if (ifs.eof() || this_frame.natm == 0) break;
			frames.push_back(this_frame);
			++frameno;
	}

	return frameno;
}

//----------------------------------------------------------------------------//

void read_xyz_frame(size_t& lineno, size_t& imcon, std::istream& is, opt::lammpstrj_frame& frame, double t, std::vector<double> box) 
{
    assert(frame.natm == 0); // empty

    const auto input_error = [&lineno](const std::string& message) {
        throw std::runtime_error("XYZ input at line " + std::to_string(lineno) +
                                 ": " + message);
    };

    std::string line;
    if (!std::getline(is, line)) {
        if (is.eof() && !is.bad())
            return; // Normal end of input between frames.
        input_error("could not read the atom-count line");
    }
    ++lineno;
    std::istringstream iss(line);
    iss >> std::ws;
    size_t atom_count = 0;
    std::string extra;
    if (iss.peek() == '-' || !(iss >> atom_count) || atom_count == 0 || iss >> extra ||
        atom_count > std::numeric_limits<size_t>::max() / (4 * sizeof(double))) {
        input_error("expected one positive integer atom count within the supported range");
    }

    ++lineno;
    if (!std::getline(is, line))
        input_error("incomplete frame: expected the XYZ comment line");

    frame.natm = atom_count;
    frame.xyz        = new double[3*frame.natm];
    frame.name       = new char[frame.natm*MAX_CHAR_LENGTH];
    frame.atomic_num = new int[frame.natm];
    frame.mass       = new double[frame.natm];
    frame.cell       = new double[9];
    frame.rcell       = new double[9];
    frame.dipind       = new double[4*frame.natm];
    frame.dipmol       = new double[4*frame.natm];
    frame.molecid       = new int[frame.natm];

    std::fill(frame.name, frame.name + frame.natm*MAX_CHAR_LENGTH, '0');
    std::fill(frame.cell, frame.cell + 9, 0);
    frame.time = t;
    
    for(size_t i = 0; i < 9; ++i){
        frame.cell[i] =  box[i];
    }

    if (imcon < 3) {

        frame.vol = frame.cell[0]*frame.cell[4]*frame.cell[8];

        std::fill(frame.rcell, frame.rcell + 9, 0.0);
        frame.rcell[0] = 1.0/frame.cell[0];
        frame.rcell[4] = 1.0/frame.cell[4];
        frame.rcell[8] = 1.0/frame.cell[8];

    } else if (imcon == 3){	// non-orthorhombic

        double cross_x = frame.cell[4]*frame.cell[8] - frame.cell[5]*frame.cell[7];
        double cross_y = -frame.cell[1]*frame.cell[8] + frame.cell[7]*frame.cell[2];
        double cross_z = frame.cell[1]*frame.cell[5] - frame.cell[4]*frame.cell[2];

        frame.vol = frame.cell[0]*cross_x + frame.cell[3]*cross_y + frame.cell[6]*cross_z; // [a] . [b] x [c]

	    calculate_rcell(lineno, frame.cell, frame.rcell);

    } else {
	std::ostringstream oss;
	oss << "Do not know how to calculate volume for imcon = " 
     	    << imcon << std::endl;
	throw std::runtime_error(oss.str());
    }



    for (size_t n = 0; n < frame.natm; ++n) {
        ++lineno;
        if (!std::getline(is, line))
            input_error("incomplete frame: expected atom " + std::to_string(n + 1) +
                        " of " + std::to_string(frame.natm));

        iss.clear();
        iss.str(line);
        std::string label;
        if (!(iss >> label))
            input_error("expected an atom label followed by three finite numeric components");
        for (size_t component = 0; component < 3; ++component) {
            std::string field;
            if (!(iss >> field))
                input_error("expected an atom label followed by three finite numeric components");
            std::istringstream number(field);
            double value;
            if (!(number >> value) || !number.eof() || !std::isfinite(value))
                input_error("invalid numeric component '" + field +
                            "'; expected a finite coordinate or velocity");
            frame.xyz[3*n + component] = value;
        }
        frame.molecid[n] = (n + 1) / 3;
        

    
    }
    
    
    iss.clear();
}


} // namespace

//----------------------------------------------------------------------------//
