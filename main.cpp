#include <iostream>
#include <vector>
#include <utility>
#include <boost/random.hpp>
#include <fstream>
#include <chrono>
#include <list>
#include <algorithm>
#include <functional>

typedef double double_t;
struct coordinate_t{
    double_t t,x,y;
    coordinate_t(){ t = 0; x = 0; y = 0; }
    coordinate_t( double_t _t, double_t _x, double_t _y ) { t = _t; x = _x; y = _y; }
};
typedef std::vector<coordinate_t> trajectory_t;
typedef std::function<coordinate_t(const coordinate_t&)> aFunc;
typedef std::function<coordinate_t(const coordinate_t&, double_t, double_t)> bFunc;

trajectory_t getTrajectory(
        double_t tBegin,
        double_t tEnd,
        aFunc a_func,
        bFunc b_func,
        coordinate_t x0,
        size_t N
        )
{
    static bool firstCall = true;
    unsigned seed = 0;
    if( firstCall )
    {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
        firstCall = false;
    }

    static boost::mt19937 generator( seed );
    static boost::normal_distribution<double_t> norm( 0., 1. );
    static boost::variate_generator< boost::mt19937, boost::normal_distribution<double> >
            normalSampler( generator, norm );

    trajectory_t answer( N, x0 );

    double_t tStep = (tEnd - tBegin) / (N - 1);

    for( size_t i = 1; i < N; ++i )
    {
        double_t curT = tStep * i + tBegin;

        double_t epsilon = normalSampler();
        double_t psi = normalSampler();

        coordinate_t flow = a_func( answer[i-1] );
        coordinate_t volat = b_func( answer[i-1], epsilon, psi);
        coordinate_t res;
        res.t = curT;
        res.x = answer[i-1].x + flow.x * tStep + volat.x * sqrt( tStep );
        res.y = answer[i-1].y + flow.y * tStep + volat.y * sqrt( tStep );
        answer[i] = res;
    }

    return answer;
}

void writeTrajextoryToStream( trajectory_t &tr, std::ostream &str )
{
    for( size_t i = 0; i < tr.size(); ++i )
        str <<
               tr[i].t <<
               ',' <<
               tr[i].x <<
               ',' << tr[i].y << std::endl;
}

int main( int argc, char** argv )
{
    if( argc != 2 )
        throw std::logic_error( "bad params" );

    const size_t NUM_OF_TRAJECTORIES = atoi( argv[1] );//20;

    aFunc a_func = [&](const coordinate_t& X) -> coordinate_t
    {
        return coordinate_t( X.t, -X.x, -X.y );
    };

    bFunc b_func = [&](const coordinate_t& X, double_t w1, double_t w2 ) -> coordinate_t
    {
        return coordinate_t( X.t, 0. w1 + 0.4 * w2, 0.4 * w1 + 0. * w2 );
    };

    std::vector<trajectory_t> trajectories;
    coordinate_t X0( 0, 1, 1 );
    for( size_t i = 0; i < NUM_OF_TRAJECTORIES; ++i )
    {
        std::fstream fStream( "res1_" + std::to_string( i ) + ".csv", std::ios::out );
        trajectory_t tr = getTrajectory(
                    0,
                    2,
                    a_func,
                    b_func,
                    X0,
                    1000
                    );
        trajectories.push_back( tr );
        writeTrajextoryToStream( tr, fStream );
    }

    trajectory_t mean( 1000 );

    for( size_t i = 0; i < 1000; ++i )
    {
        mean[i].t = 0.002 * i;
        mean[i].x = 0;
        mean[i].y = 0;
        for( size_t j = 0; j < trajectories.size(); ++j )
        {
            mean[i].x += trajectories[j][i].x / trajectories.size();
            mean[i].y += trajectories[j][i].y / trajectories.size();
        }
    }

    std::fstream meanStr( "mean1.csv", std::ios::out );

    writeTrajextoryToStream( mean, meanStr );

    return EXIT_SUCCESS;
}
