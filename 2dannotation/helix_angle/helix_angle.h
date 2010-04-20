#ifndef _helix_angle_h_
#define _helix_angle_h_

class HelixAngle
{
public:
	// LIFECYCLE
	HelixAngle(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;
private:
	void readOptions(int argc, char* argv[]);

	mccore::Molecule* loadFile (const std::string &astrFilename);
};

#endif /*_helix_angle_h_*/
