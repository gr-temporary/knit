
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string.h>
#include "CImg.h"

using namespace cimg_library;

typedef int Chromosome;
typedef std::vector<Chromosome> DNA;
typedef std::pair<int, int> Point;
typedef CImg<float> Image;

double uniformRandom()
{
	return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
}

int intRandom(int start, int end) {
	return start + rand() % (end - start);
}
// return a normally distributed random number
double normalRandom()
{
	double u1 = uniformRandom();
	double u2 = uniformRandom();
	return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

bool inCircle(int x, int y, int r) {
	int dx = x - r;
	int dy = y - r;
	return (dx * 2 + 1) * (dx * 2 + 1) + (dy * 2 + 1) * (dy * 2 + 1) <= (r * 2 + 1) * (r * 2 + 1);
}

float normalizeImage(Image &canvas) {
	float mean = 0.0;
	int c = 0;
	int r = canvas.width() / 2;
	for (int i = 0; i < canvas.width(); i++) {
		for (int j = 0; j < canvas.height(); j++) {
			if (inCircle(i, j, r)) {
				mean += canvas(i, j, 0, 0);
				c++;
			}
		}
	}
	mean /= c;
	float dev = 0.0;
	for (int i = 0; i < canvas.width(); i++) {
		for (int j = 0; j < canvas.height(); j++) {
			if (inCircle(i, j, r)) {
				float d = canvas(i, j, 0, 0) - mean;
				dev += pow(d, 2);
				canvas(i, j, 0, 0) = d;
			}
		}
	}
	dev = sqrt(dev);
	for (int i = 0; i < canvas.width(); i++) {
		for (int j = 0; j < canvas.height(); j++) {
			if (inCircle(i, j, r)) {
				canvas(i, j, 0, 0) /= dev;
			}
		}
	}
	for (int i = 0; i < canvas.width(); i++) {
		for (int j = 0; j < canvas.height(); j++) {
			if (!inCircle(i, j, r)) {
				canvas(i, j, 0, 0) = 0.0;
			}
		}
	}
	return dev;
}

struct Genome {
	DNA dna;
	double fitness;
	double weight;

	void init(int size, int nails) {
		dna.resize(size);
		int minOffset = nails * 1 / 5;
		int maxOffset = nails * 4 / 5;
		dna[0] = rand() % nails;
		for (int i = 1; i < size; i++) {
			dna[i] = (dna[i - 1] + intRandom(minOffset, maxOffset)) % nails;
		}
		fitness = 0.0;
	}

	void mutate(double spread, int nails) {
		for (int i = 0; i < dna.size(); i++) {
			int offset = floor(pow(normalRandom() * spread, 2) + 0.5);
			offset += nails; // to eliminate negative offsets and %
			dna[i] = (dna[i] + offset) % nails;
		}
	}

	void draw(Image &canvas, const std::vector<Point> &nails, const float threadOpacity, const std::vector<Point> &bounds) {
		canvas.fill(0.0f);
		for (int i = 0; i < dna.size() - 1; i++) {
			int x0 = nails[dna[i]].first;
			int x1 = nails[dna[i + 1]].first;
			int y0 = nails[dna[i]].second;
			int y1 = nails[dna[i + 1]].second;
			int dx = cimg::abs(x1 - x0);
			int dy = cimg::abs(y1 - y0);
			int sx = (x0 < x1) ? 1 : -1;
			int sy = (y0 < y1) ? 1 : -1;
			int err = dx - dy;
			while (true) {
				float c = canvas(x0, y0, 0, 0);
				c += threadOpacity;
				if (c > 1.0f) c = 1.0f;
				canvas(x0, y0, 0, 0) = c;				

				if ((x0 == x1) && (y0 == y1)) break;
				int e2 = 2 * err;
				if (e2 >-dy) { err -= dy; x0 += sx; }
				if (e2 < dx) { err += dx; y0 += sy; }
			}
		}
		normalizeImage(canvas);
	}

	// ORDER BY fitness DESC
	bool operator <(const Genome &g) {
		return fitness > g.fitness;
	}
};

typedef std::vector<Genome> Genomes;

double convolve(const Image &kernel, const Image &canvas) {
	double result = 0.0;
	double n = 0;
	int r = canvas.width() / 2;
	for (int i = 0; i < canvas.width(); i++) {
		for (int j = 0; j < canvas.height(); j++) {
			if (inCircle(i, j, r)) {
				result += kernel(i, j, 0, 0) * canvas(i, j, 0, 0);
				n++;
			}
		}
	}
	return result;
}

struct Population {
	Genomes population;
	double weightSum;

	void init(int populationSize, int genomeSize, int nails) {
		population.resize(populationSize);
		for (int i = 0; i < population.size(); i++) {
			population[i].init(genomeSize, nails);
		}
	}

	double calculateFitness(const Image &kernel, Image &canvas, const std::vector<Point> &nails, const float threadOpacity, const std::vector<Point> &bounds) {
		double minFitness = 10e+7;
		double maxFitness = -10e+7;
		double average = 0.0;

		for (int i = 0; i < population.size(); i++) {
			//printf("%i\n", i);
			population[i].draw(canvas, nails, threadOpacity, bounds);
			double fitness = convolve(kernel, canvas);
			population[i].fitness = fitness;
			if (fitness > maxFitness) maxFitness = fitness;
			if (fitness < minFitness) minFitness = fitness;
			average += fitness;
		}
		std::sort(population.begin(), population.end());

		weightSum = 0.0;
		for (int i = 0; i < population.size(); i++) {
			double weight = (population[i].fitness - minFitness) / (maxFitness - minFitness);
			population[i].weight = weightSum;
			weightSum += weight;
		}
		
		return average / population.size();
	}

	Genome* getParent() {
		double random = uniformRandom() * weightSum / 2; // only first (successful) part, because fuck you
		int i;
		for (i = population.size() - 1; i > 0; i--) {
			if (random > population[i].weight)
				break;
		}
		//printf("Parent: %i (%lf)\n", i, population[i].fitness);
		return &(population[i]);
	}

	Genome* getParentOrdinal() {
		double random = pow(uniformRandom(), 2) * population.size() / 2;
		int i = floor(random);
		//printf("%i ", i);
		return &(population[i]);
	}
};

int getCrossoverJumpSize(int size) {
	return floor(uniformRandom() * size / 2 + size / 3);
}

Genome crossover(Genome *mother, Genome *father) {
	Genome child;
	int size = mother->dna.size();

	child.dna.resize(size);

	const Genome *donor;
	int i = 0;
	int j = 0;
	int jump = getCrossoverJumpSize(size);
	donor = uniformRandom() < 0.5 ? mother : father;
	while (i < size) {
		child.dna[i] = donor->dna[i];
		i++;
		j++;
		if (j == jump) {
			j = 0;
			jump = getCrossoverJumpSize(size);
			donor = donor == father ? mother : father;
		}
	}
	return child;
}

void generateNewPopulation(Population &base, Population &result, int nails) {
	for (int i = 0; i < result.population.size(); i++) {
		// select 2 parents
		Genome* mother = base.getParentOrdinal();
		Genome* father = base.getParentOrdinal();

		// crossover
		result.population[i] = crossover(mother, father);

		// mutate
		result.population[i].mutate(1.5, nails);
	}
}

Image loadImage(char *filename, int imageSize) {
	Image result(imageSize, imageSize, 1, 1, 0.0);
	CImg<unsigned char> image(filename);

	image.resize(imageSize, imageSize);

	for (int i = 0; i < imageSize; i++) {
		for (int j = 0; j < imageSize; j++) {
			double c = image(i, j, 0, 0);
			c = 1.0 - c / 255.0;
			result(i, j, 0, 0) = c;

			// can't harm
			if (!inCircle(i, j, imageSize / 2)) {
				result(i, j, 0, 0) = 0.0;
			}
		}
	}

	return result;
}

void drawHistory(std::vector<double> history) {
	Image img(history.size(), 200, 1, 1, 0.0);
	std::vector<double>::iterator max = std::max_element(history.begin(), history.end());
	std::vector<double>::iterator min = std::min_element(history.begin(), history.end());
	
}

void run() {
	int nails = 150;
	int populationSize = 100;
	int genomeSize = nails * 4;
	int imageSize = 100;
	int iterations = 50;
	int stableIterations = 5;
	float threadOpacity = 1.0 / 10;
	double threshold = 1e-7;
	char filename[256] = "test-i.png";

	FILE *config = NULL;
	fopen_s(&config, "config.txt", "r");
	
	char line[256];
	while (fscanf_s(config, "%s", line, 256) != EOF) {
		if (strcmp(line, "file") == 0) {
			//strcpy_s(filename, line);
			fscanf_s(config, "%s", filename, 256);
		} else
		if (strcmp(line, "nails") == 0) {
			fscanf_s(config, "%i", &nails);
		} else
		if (strcmp(line, "populationsize") == 0) {
			fscanf_s(config, "%i", &populationSize);
		} else
		if (strcmp(line, "genomesize") == 0) {
			fscanf_s(config, "%i", &genomeSize);
		} else
		if (strcmp(line, "imagesize") == 0) {
			fscanf_s(config, "%i", &imageSize);
		} else
		if (strcmp(line, "iterations") == 0) {
			fscanf_s(config, "%i", &iterations);
		}
		if (strcmp(line, "stableiterations") == 0) {
			fscanf_s(config, "%i", &stableIterations);
		}
		else
		if (strcmp(line, "threshold") == 0) {
			fscanf_s(config, "%f", &threshold);
		}
		else
		if (strcmp(line, "opacity") == 0) {
			fscanf_s(config, "%f", &threadOpacity);
		} else {
			
		}
	}

	fclose(config);

	std::vector<Point> nailPoints;
	nailPoints.resize(nails);
	for (int i = 0; i < nails; i++) {
		double a = 3.14159265359 * 2 * i / nails;
		double s = std::sin(a);
		double c = std::cos(a);
		double r = imageSize / 2 - 1;
		nailPoints[i].first = floor(imageSize / 2 + c * r);
		nailPoints[i].second = floor(imageSize / 2 + s * r);
	}
	std::vector<Point> bounds;
	bounds.resize(imageSize);
	for (int i = 0; i < imageSize; i++) {
		double hw = imageSize / 2 - 0.5;
		int off = sqrt(hw * hw - pow(i - hw, 2));
		bounds[i].first = floor(imageSize / 2 - off);
		bounds[i].second = floor(imageSize / 2 + off);
	}

	Population populations[2];
	populations[0].init(populationSize, genomeSize, nails);
	populations[1].init(populationSize, genomeSize, nails);

	Image kernel = loadImage(filename, imageSize);
	Image canvas(imageSize, imageSize, 1, 1, 0.0f);

	normalizeImage(kernel);

	double fitness = 0.0;
	double prevFitness = 1.0;
	int stableGeneration = 0;
	int generation = 0;
	bool enough = false;
	std::vector<double> history;
	history.reserve(iterations);
	
	while (!enough && generation < iterations) {
		int current = generation % 2;
		int next = (generation + 1) % 2;

		prevFitness = fitness;
		fitness = populations[current].calculateFitness(kernel, canvas, nailPoints, threadOpacity, bounds);

		history.push_back(fitness);
		printf("Generation: %i,\tFitness: %.7lf,\tMax: %.7lf\tMin: %.7lf\n", generation, fitness, populations[current].population.front().fitness, populations[current].population.back().fitness);

		generateNewPopulation(populations[current], populations[next], nails);

		if (abs(1.0 - fitness / prevFitness) < threshold) {
			stableGeneration++;
		}
		else {
			stableGeneration = 0;
		}
		if (stableGeneration >= stableIterations) {
			enough = true;
		}
		generation++;
	}
	fitness = populations[generation % 2].calculateFitness(kernel, canvas, nailPoints, threadOpacity, bounds);
	printf("Generation: %i,\tFitness: %lf\n", generation, fitness);

	Genome *win = &(populations[generation % 2].population.front());
	win->draw(canvas, nailPoints, threadOpacity, bounds);
	canvas.normalize(0, 255);
	canvas *= -1;
	canvas += 255;
	canvas.display("Winner");

	win = &(populations[generation % 2].population.back());
	win->draw(canvas, nailPoints, threadOpacity, bounds);
	canvas.normalize(0, 255);
	canvas *= -1;
	canvas += 255;
	canvas.display("Not winner");
}

int main(int argc, char **argv) {
	srand(clock());
	/*for (int i = 0; i < 20; i++) {
		printf("%i\n", (int)floor(normalRandom() * 2 + 0.5));
	}*/

	/*Image test1 = loadImage("stripes.png", 100);
	Image test2(100, 100, 1, 1, 0.0);
	for (int i = 0; i < 100 - 3; i++) {
		for (int j = 0; j < 100; j++) {
			test2(i, j, 0, 0) = test1(i + 3, j, 0, 0);
		}
	}
	normalizeImage(test1);
	normalizeImage(test2);
	test1.display();
	test2.display();
	double r = convolve(test1, test2);
	printf("%lf", r);
	test1 -= test2;
	test1.display();*/

	run();

	char c;
	//scanf_s("%c", &c);
}