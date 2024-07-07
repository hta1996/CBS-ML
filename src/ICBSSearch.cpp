#include "ICBSSearch.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "RectangleReasoning.h"

#define MAX_RUNTIME4PAIR 6000
// takes the paths_found_initially and UPDATE all (constrained) paths found for agents from curr to start
// also, do the same for ll_min_f_vals and paths_costs (since its already "on the way").
inline void ICBSSearch::updatePaths(ICBSNode* curr)
{
    for(int i = 0; i < num_of_agents; i++)
        paths[i] = &paths_found_initially[i];
    vector<bool> updated(num_of_agents, false);  // initialized for false

    while (curr->parent != NULL)
    {
        if (!updated[curr->agent_id])
        {
            if (curr->path.empty())
                cout << "ERROR" << endl;
            paths[curr->agent_id] = &(curr->path);
            updated[curr->agent_id] = true;
        }
        curr = curr->parent;
    }
}

// return the minimal length of the path
int ICBSSearch::collectConstraints(ICBSNode* curr, int agent_id, std::vector <std::list< std::pair<int, int> > >& cons_vec)
{
    std::clock_t t1 = std::clock();
    // extract all constraints on agent_id
    int minLength = 0;
    list < tuple<int, int, int> > constraints;
    int max_timestep = -1;
    while (curr != dummy_start)
    {
        if (curr->agent_id == agent_id)
        {
            for (Constraint constraint : curr->constraints)
            {
                constraints.push_back(constraint);
                if (get<2>(constraint) > max_timestep) // calc constraints' max_timestep
                    max_timestep = get<2>(constraint);
                if(get<2>(constraint) > minLength && get<0>(constraint) < 0) // this is a path length constraint
                    minLength = get<2>(constraint);
                else if (get<2>(constraint) >= minLength &&
                    get<0>(constraint) == paths[agent_id]->back().location && get<1>(constraint) < 0) // or this is a vertex constraint at the goal location
                    minLength = get<2>(constraint) + 1;
            }
        }
        curr = curr->parent;
    }
    for (Constraint constraint : initial_constraints[agent_id])
    {
        constraints.push_back(constraint);
        if (get<2>(constraint) > max_timestep) // calc constraints' max_timestep
            max_timestep = get<2>(constraint);
        if (get<2>(constraint) > minLength && get<0>(constraint) < 0) // this is a path length constraint
            minLength = get<2>(constraint);
        else if (get<2>(constraint) >= minLength &&
            get<0>(constraint) == paths[agent_id]->back().location && get<1>(constraint) < 0) // or this is a vertex constraint at the goal location
            minLength = get<2>(constraint) + 1;
    }

    // initialize a constraint vector of length max_timestep+1. Each entry is an empty list< pair<int,int> > (loc1,loc2)
    //std::vector <std::list< std::pair<int, int> > > cons_vec;
    cons_vec.resize(max_timestep + 1);
    
    for (list<Constraint>::iterator it = constraints.begin(); it != constraints.end(); it++)
    {
        if (get<0>(*it) < 0) // barrier constraint
        {
            int x1 = (-get<0>(*it) - 1) / ml->cols, y1 = (-get<0>(*it) - 1) % ml->cols;
            int x2 = get<1>(*it) / ml->cols, y2 = get<1>(*it) % ml->cols;
            if (x1 == x2)
            {
                if (y1 < y2)
                    for (int i = 0; i <= y2 - y1; i++)
                        cons_vec[get<2>(*it) - i].push_back(make_pair(x1 *  ml->cols + y2 - i, -1));
                else
                    for (int i = 0; i <= y1 - y2; i++)
                        cons_vec[get<2>(*it) - i].push_back(make_pair(x1 *  ml->cols + y2 + i, -1));
            }
            else // y1== y2
            {
                if (x1 < x2)
                    for (int i = 0; i <= x2 - x1; i++)
                        cons_vec[get<2>(*it) - i].push_back(make_pair((x2 - i) *  ml->cols + y1, -1));
                else
                    for (int i = 0; i <= x1 - x2; i++)
                        cons_vec[get<2>(*it) - i].push_back(make_pair((x2 + i) *  ml->cols + y1, -1));
            }
        }
        else
            (cons_vec)[get<2>(*it)].push_back(make_pair(get<0>(*it), get<1>(*it)));
    }
    runtime_updatecons += (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;
    return minLength;

}



int ICBSSearch::computeHeuristics(ICBSNode& curr)
{
    // create conflict graph
    vector<int> CG(num_of_agents * num_of_agents, 0);
    int num_of_CGedges = 0;
    if (h_type == heuristics_type::DG || h_type == heuristics_type::WDG)
    {
        bool succeed = buildDependenceGraph(curr);
        if (!succeed)
            return -1;
        for (int i = 0; i < num_of_agents; i++)
        {
            for (int j = i + 1; j < num_of_agents; j++)
            {
                auto got = curr.conflictGraph.find(i * num_of_agents + j);
                if (got != curr.conflictGraph.end() && got->second > 0)
                {
                    CG[i * num_of_agents + j] = got->second;
                    CG[j * num_of_agents + i] = got->second;
                    num_of_CGedges++;
                }
            }
        }
    }
    else
    {
        for (list<std::shared_ptr<Conflict>>::const_iterator it = curr.cardinalConf.begin(); it != curr.cardinalConf.end(); ++it)
        {
            if (!CG[get<0>(**it) * num_of_agents + get<1>(**it)])
            {
                CG[get<0>(**it) * num_of_agents + get<1>(**it)] = true;
                CG[get<1>(**it) * num_of_agents + get<0>(**it)] = true;
                num_of_CGedges++;
            }
        }
    }

    clock_t t = std::clock();
    int rst;
    if (h_type == heuristics_type::WDG)
    {
        rst = weightedVertexCover(CG, num_of_agents);
    }
    else
    {
        // Minimum Vertex Cover
        if (curr.parent == NULL) // root node of CBS tree
            rst = minimumVertexCover(CG, -1, num_of_agents, num_of_CGedges);
        else
            rst = minimumVertexCover(CG, curr.parent->h_val, num_of_agents, num_of_CGedges);
    }
    mvc_runtime += (std::clock() - t) * 1000.0 / CLOCKS_PER_SEC;

    return rst;
}

int ICBSSearch::getEdgeWeight(int a1, int a2, const vector<list<Constraint>> & constraints, ICBSNode& node, bool cardinal, bool& hit)
{
    HTableEntry newEntry(a1, a2, &node);
    if (h_type != heuristics_type::CG)
    {
        time_t t = std::clock();
        HTable::const_iterator got = hTable[a1][a2].find(newEntry);
        bookingSearchtime += (std::clock() - t) * 1000.0/ CLOCKS_PER_SEC;

        if (got != hTable[a1][a2].end())
        {
            bookingHitTimes++;
            hit = true;
            return got->second;
        }
        
    }
    hit = false;
    vector<list<Constraint>> cons(2);
    cons[0] = constraints[a1];
    cons[1] = constraints[a2];
    int cost_shortestPath = paths[a1]->size() + paths[a2]->size() - 2;
    runtime = (std::clock() - start) * 1000.0 / CLOCKS_PER_SEC;
    int scr = 0;
    if(screen == 2)
    {
        scr = 1;
        std::cout << "Agents " << a1 << " and " << a2 << " in node " << node.time_generated << " : ";
    }
    int rst;
    if (cardinal)
        rst = 1;
    else if (h_type == heuristics_type::DG || h_type == heuristics_type::WDG)
    {
        // get mdds
        
        const MDD* mdd1 = buildMDD(node, a1);
        const MDD* mdd2 = buildMDD(node, a2);

        /*ICBSNode* curr = &node;
        while (curr != dummy_start)
        {
            if (mdd1 == NULL && curr->agent_id == a1)
            {
                if (curr->mdd == NULL)
                {
                    buildMDD(*curr, a1);
                }
                mdd1 = curr->mdd;
            }
            else if (mdd2 == NULL && curr->agent_id == a2)
            {
                if (curr->mdd == NULL)
                {
                    buildMDD(*curr, a2);
                }
                mdd2 = curr->mdd;
            }
            curr = curr->parent;
        }
        if (mdd1 == NULL)
        {
            if (mdds_initially[a1] == NULL)
            {
                buildMDD(*dummy_start, a1);
            }
            mdd1 = mdds_initially[a1];
        }
        if (mdd2 == NULL)
        {
            if (mdds_initially[a2] == NULL)
            {
                buildMDD(*dummy_start, a2);
            }
            mdd2 = mdds_initially[a2];
        }*/
        if (mdd1->levels.size() > mdd2->levels.size()) // swap
        {
            const MDD* temp = mdd1;
            mdd1 = mdd2;
            mdd2 = temp;
        }
        if(!SyncMDDs(*mdd1, *mdd2))
            rst = 1;
        else
            rst = 0;
        if (screen == 2)
        {
            std::cout << rst << std::endl;
            //test
            vector<SingleAgentICBS*> engines(2);
            engines[0] = search_engines[a1];
            engines[1] = search_engines[a2];
            vector<vector<PathEntry>> initial_paths(2);
            initial_paths[0] = *paths[a1];
            initial_paths[1] = *paths[a2];
            ICBSSearch solver2(ml, engines, cons, initial_paths, 1.0, 0, heuristics_type::CG, true, rectangleReasoning, INT_MAX, std::min(1000.0, time_limit - runtime), 0);
            solver2.runICBSSearch();
            if ((rst ==0 && solver2.solution_cost > (int)(initial_paths.size() - initial_paths.size() + 2)) &&
                (rst == 1 && solver2.solution_cost == (int)(initial_paths.size() - initial_paths.size() + 2)))
            {
                cout << "SyncMDD WRONG" << endl;
                exit(10);
            }
        }
    }
    if (h_type == heuristics_type::WDG && rst > 0)
    {
        vector<SingleAgentICBS*> engines(2);
        engines[0] = search_engines[a1];
        engines[1] = search_engines[a2];
        vector<vector<PathEntry>> initial_paths(2);
        initial_paths[0] = *paths[a1];
        initial_paths[1] = *paths[a2];
        double cutoffTime = std::min(MAX_RUNTIME4PAIR * 1.0, time_limit - runtime);
        int upperbound = initial_paths[0].size() + initial_paths[1].size() + 10;
        ICBSSearch solver(ml, engines, cons, initial_paths, 1.0, max(rst, 0), heuristics_type::CG, true, rectangleReasoning, upperbound, cutoffTime, scr);
        solver.max_num_of_mdds = this->max_num_of_mdds;
        solver.runICBSSearch();
        if (solver.runtime >= cutoffTime) // time out
            rst = solver.min_f_val - cost_shortestPath; // using lowverbound to approximate
        else if (solver.solution_cost  < 0) // no solution
            rst = solver.solution_cost;
        else
            rst = solver.solution_cost - cost_shortestPath;
    }
    hTable[a1][a2][newEntry] = rst;
    return rst;
}

bool ICBSSearch::buildDependenceGraph(ICBSNode& node)
{
    // extract all constraints
    vector<list<Constraint>> constraints = initial_constraints;
    ICBSNode* curr = &node;
    while (curr != dummy_start)
    {
        for (auto constraint : curr->constraints)
        {
            constraints[curr->agent_id].push_back(constraint);
        }
        curr = curr->parent;
    }
    
    if (screen == 2)
    {
        for (size_t i = 0; i < constraints.size(); i++)
        {
            if (constraints[i].empty())
                continue;
            std::cout << "Constraints for agent " << i << ":";
            for (auto constraint: constraints[i])
                std::cout << constraint;
            std::cout <<std::endl;
        }
    }

    for (auto conflict : node.cardinalConf)
    {
        int a1 = min(get<0>(*conflict), get<1>(*conflict));
        int a2 = max(get<0>(*conflict), get<1>(*conflict));
        int idx = a1 * num_of_agents + a2;
        if (h_type == heuristics_type::DG)
        {
            node.conflictGraph[idx] = 1;
        }
        else if (node.conflictGraph.find(idx) == node.conflictGraph.end())
        {
            time_t t = std::clock();
            bool hit;
            int w = getEdgeWeight(a1, a2, constraints, node, true, hit);
            cardinal_pair_runtime += (std::clock() - t) * 1000.0 / CLOCKS_PER_SEC;
            if(!hit)
                cardinal_pair_num ++;
            if (w < 0) // no solution
                return false;
            
            node.conflictGraph[idx] = w;
            if (w > 1 && !hit)
            {
                cardinal_heuristic_num++;
                cardinal_heuristic_value += w;
            }
        }
    }
    for (auto conflict : node.semiConf)
    {
        int a1 = min(get<0>(*conflict), get<1>(*conflict));
        int a2 = max(get<0>(*conflict), get<1>(*conflict));
        int idx = a1 * num_of_agents + a2;
        if (node.conflictGraph.find(idx) == node.conflictGraph.end())
        {
            time_t t = std::clock();
            bool hit;
            int w = getEdgeWeight(a1, a2, constraints, node, false, hit);
            not_cardinal_pair_runtime += (std::clock() - t) * 1000.0 / CLOCKS_PER_SEC;
            if (!hit)
                not_cardinal_pair_num++;
            if (w < 0) //no solution
                return false;
            node.conflictGraph[idx] = w;
            if (w > 0 && !hit)
            {
                not_cardinal_heuristic_num++;
                not_cardinal_heuristic_value += w;
            }
        }
    }
    for (auto conflict : node.nonConf)
    {
        int a1 = min(get<0>(*conflict), get<1>(*conflict));
        int a2 = max(get<0>(*conflict), get<1>(*conflict));
        int idx = a1 * num_of_agents + a2;
        if (node.conflictGraph.find(idx) == node.conflictGraph.end())
        {
            time_t t = std::clock();
            bool hit;
            int w = getEdgeWeight(a1, a2, constraints, node, false, hit);
            not_cardinal_pair_runtime += (std::clock() - t) * 1000.0 / CLOCKS_PER_SEC;
            if (!hit)
                not_cardinal_pair_num++;
            if (w < 0) // no solution
                return false;
            node.conflictGraph[idx] = w;
            if (w > 0)
            {
                not_cardinal_heuristic_num++;
                not_cardinal_heuristic_value += w;
            }
        }
    }
    return true;
}


// deep copy of all conflicts except ones that involve the particular agent
// used for copying conflicts from the parent node to the child nodes
void ICBSSearch::copyConflicts(const std::list<std::shared_ptr<Conflict>>& conflicts,
    std::list<std::shared_ptr<Conflict>>& copy, int excluded_agent) const
{
    for (std::list<std::shared_ptr<Conflict>>::const_iterator it = conflicts.begin(); it != conflicts.end(); ++it)
    {
        if (get<0>(**it) != excluded_agent && get<1>(**it) != excluded_agent)
        {
            copy.push_back(*it);
        }
    }
}


void ICBSSearch::findConflicts(ICBSNode& curr)
{
    if (curr.parent != NULL)
    {
        // Copy from parent
        copyConflicts(curr.parent->rectSemiConf, curr.rectSemiConf, curr.agent_id);
        copyConflicts(curr.parent->rectNonConf, curr.rectNonConf, curr.agent_id);
        copyConflicts(curr.parent->cardinalConf, curr.cardinalConf, curr.agent_id);
        copyConflicts(curr.parent->semiConf, curr.semiConf, curr.agent_id);
        copyConflicts(curr.parent->nonConf, curr.nonConf, curr.agent_id);
        copyConflicts(curr.parent->unknownConf, curr.unknownConf, curr.agent_id);

        // detect new conflicts
        int a1 = curr.agent_id;
        for (int a2 = 0; a2 < num_of_agents; a2++)
        {
            if(a1 == a2)
                continue;
            else if (search_engines[a1]->num_of_conf == 0) // New path does not have conflicts with others before it reaches its goal
            {
                if (paths[a1]->size() + 1 < paths[a2]->size())
                {
                    int loc1 = paths[a1]->back().location;
                    for (size_t timestep = paths[a1]->size(); timestep < paths[a2]->size(); timestep++)
                    {
                        int loc2 = paths[a2]->at(timestep).location;
                        if (loc1 == loc2)
                        {
                            curr.unknownConf.push_front(std::shared_ptr<Conflict>(new Conflict(a1, a2, loc1, -1, timestep))); // It's at least a semi conflict
                        }
                    }
                }
                continue;
            }
            size_t min_path_length = paths[a1]->size() < paths[a2]->size() ? paths[a1]->size() : paths[a2]->size();
            for (size_t timestep = 0; timestep < min_path_length; timestep++)
            {
                int loc1 = paths[a1]->at(timestep).location;
                int loc2 = paths[a2]->at(timestep).location;
                if (loc1 == loc2)
                {
                    curr.unknownConf.push_back(std::shared_ptr<Conflict>(new Conflict(a1, a2, loc1, -1, timestep)));
                }
                else if (timestep < min_path_length - 1
                    && loc1 == paths[a2]->at(timestep + 1).location
                    && loc2 == paths[a1]->at(timestep + 1).location)
                {
                    curr.unknownConf.push_back(std::shared_ptr<Conflict>(new Conflict(a1, a2, loc1, loc2, timestep + 1))); // edge conflict
                }
            }
            if (paths[a1]->size() != paths[a2]->size())
            {
                int a1_ = paths[a1]->size() < paths[a2]->size() ? a1 : a2;
                int a2_ = paths[a1]->size() < paths[a2]->size() ? a2 : a1;
                int loc1 = paths[a1_]->back().location;
                for (size_t timestep = min_path_length; timestep < paths[a2_]->size(); timestep++)
                {
                    int loc2 = paths[a2_]->at(timestep).location;
                    if (loc1 == loc2)
                    {
                        curr.unknownConf.push_front(std::shared_ptr<Conflict>(new Conflict(a2_, a1_, loc1, -1, timestep))); // It's at least a semi conflict
                    }
                }
            }
        }
    }
    else
    {
        for(int a1 = 0; a1 < num_of_agents ; a1++)
        {
            for (int a2 = a1 + 1; a2 < num_of_agents; a2++)
            {
                size_t min_path_length = paths[a1]->size() < paths[a2]->size() ? paths[a1]->size() : paths[a2]->size();
                for (size_t timestep = 0; timestep < min_path_length; timestep++)
                {
                    int loc1 = paths[a1]->at(timestep).location;
                    int loc2 = paths[a2]->at(timestep).location;
                    if (loc1 == loc2)
                    {
                        curr.unknownConf.push_back(std::shared_ptr<Conflict>(new Conflict(a1, a2, loc1, -1, timestep)));
                    }
                    else if (timestep < min_path_length - 1
                        && loc1 == paths[a2]->at(timestep + 1).location
                        && loc2 == paths[a1]->at(timestep + 1).location)
                    {
                        curr.unknownConf.push_back(std::shared_ptr<Conflict>(new Conflict(a1, a2, loc1, loc2, timestep + 1)));
                    }
                }
                if (paths[a1]->size() != paths[a2]->size())
                {
                    int a1_ = paths[a1]->size() < paths[a2]->size() ? a1 : a2;
                    int a2_ = paths[a1]->size() < paths[a2]->size() ? a2 : a1;
                    int loc1 = paths[a1_]->back().location;
                    for (size_t timestep = min_path_length; timestep < paths[a2_]->size(); timestep++)
                    {
                        int loc2 = paths[a2_]->at(timestep).location;
                        if (loc1 == loc2)
                        {
                            curr.unknownConf.push_front(std::shared_ptr<Conflict>(new Conflict(a2_, a1_, loc1, -1, timestep))); // It's at least a semi conflict
                        }
                    }
                }
            }
        }
    }
    curr.num_of_collisions = curr.unknownConf.size() + curr.cardinalConf.size() + curr.rectSemiConf.size() + curr.semiConf.size() + curr.rectNonConf.size() + curr.nonConf.size();
}

MDD * ICBSSearch::buildMDD(ICBSNode& node, int id)
{
    if(!mddTable.empty())
    {
        ConstraintsHasher c(id, &node);
        MDDTable::const_iterator got = mddTable[c.a].find(c);
        if (got != mddTable[c.a].end())
        {
            return got->second;
        }
    }
    MDD * mdd = new MDD();

    vector < list< pair<int, int> > > cons_vec;
    collectConstraints(&node, id, cons_vec);
    time_t t = std::clock();
    mdd->buildMDD(cons_vec,  paths[id]->size(), *search_engines[id]);
    build_mdds_runtime += (std::clock() - t) * 1000.0 / CLOCKS_PER_SEC;
    if (!mddTable.empty())
    {
        ConstraintsHasher c(id, &node);
        mddTable[c.a][c] = mdd;
    }
    return mdd;
}


void ICBSSearch::chooseConflict(ICBSNode &node)
{
    clock_t t1=std::clock();
    if (screen == 3)
        printConflicts(node);
    
    if (!node.cardinalConf.empty()) // choose the largest-cost cardinal. Break ties by the earliest cardinal
    {
        int maxWeight = -1;
        int minTime = INT_MAX;
        for (std::shared_ptr<Conflict> conflict: node.cardinalConf)
        {
            int weight = 2;
            if (get<4>(*conflict) >= (int)paths[get<0>(*conflict)]->size())
                weight = get<4>(*conflict) - paths[get<0>(*conflict)]->size() + 3;
            else if (get<4>(*conflict) >= (int)paths[get<1>(*conflict)]->size())
                weight = get<4>(*conflict) - paths[get<1>(*conflict)]->size() + 3;
            int time;
            if (get<2>(*conflict) < 0) // rectangle conflict
                time = getRectangleTime(*conflict, paths, ml->cols);
            else // vertex/edge conflict
                time = std::get<4>(*conflict);
            if (maxWeight < weight || (maxWeight == weight && time < minTime))
            {
                node.conflict = conflict;
                maxWeight = weight;
                minTime = time;
            }
        }
    }
    else if (!node.rectSemiConf.empty()) // Choose the earliest semi rect
    {
        int minTime = INT_MAX;
        for (std::shared_ptr<Conflict> conflict: node.rectSemiConf)
        {
            int time = getRectangleTime(*conflict, paths, ml->cols);
            if (time < minTime)
            {
                node.conflict = conflict;
                minTime = time;
            }
        }
    }
    else if (!node.semiConf.empty()) // Choose the earliest semi
    {
        node.conflict = node.semiConf.front();
        for (std::shared_ptr<Conflict> conflict : node.semiConf)
            if (std::get<4>(*conflict) < std::get<4>(*node.conflict))
                node.conflict = conflict;
    }
    else if (!node.rectNonConf.empty()) // Choose the earliest non rect
    {
        int minTime = INT_MAX;
        for (std::shared_ptr<Conflict> conflict : node.rectNonConf)
        {
            int time = getRectangleTime(*conflict, paths, ml->cols);
            if (time < minTime)
            {
                node.conflict = conflict;
                minTime = time;
            }
        }
    }
    else if(!node.nonConf.empty())// Choose the earliest non
    {
        node.conflict = node.nonConf.front();
        for (std::shared_ptr<Conflict> conflict : node.nonConf)
            if (std::get<4>(*conflict) < std::get<4>(*node.conflict))
                node.conflict = conflict;
    }
    else if (!node.unknownConf.empty()) // choose the earliest
    {
        node.conflict = node.unknownConf.front();
        for (std::shared_ptr<Conflict> conflict : node.unknownConf)
            if (std::get<4>(*conflict) < std::get<4>(*node.conflict))
                node.conflict = conflict;
    }
    chooseConflictTime+=(std::clock() - t1) * 1000.0/ CLOCKS_PER_SEC;
}

void ICBSSearch::classifyConflicts(ICBSNode &parent)
{
    // Classify all conflicts in unknownConf
    while (!parent.unknownConf.empty())
    {
        std::shared_ptr<Conflict> con = parent.unknownConf.front();
        parent.unknownConf.pop_front();

        bool cardinal1 = false, cardinal2 = false;
        if (get<4>(*con) >= (int)paths[get<0>(*con)]->size())
            cardinal1 = true;
        else if (paths[get<0>(*con)]->at(0).width!=1)
        {
            MDD* mdd = buildMDD(parent, get<0>(*con));
            for (size_t i = 0; i < mdd->levels.size(); i++)
                paths[get<0>(*con)]->at(i).width = mdd->levels[i].size();
            if (mddTable.empty())
                delete mdd;
        }
        if (get<4>(*con) >= (int)paths[get<1>(*con)]->size())
            cardinal2 = true;
        else if (paths[get<1>(*con)]->at(0).width!=1)
        {
            MDD* mdd = buildMDD(parent, get<1>(*con));
            for (size_t i = 0; i < mdd->levels.size(); i++)
                paths[get<1>(*con)]->at(i).width = mdd->levels[i].size();
            if (mddTable.empty())
                delete mdd;
        }

        if (get<3>(*con) >= 0) // Edge conflict
        {
            cardinal1 = paths[get<0>(*con)]->at(get<4>(*con)).width==1 && paths[get<0>(*con)]->at(get<4>(*con) - 1).width==1;
            cardinal2 = paths[get<1>(*con)]->at(get<4>(*con)).width==1 && paths[get<1>(*con)]->at(get<4>(*con) - 1).width==1;
        }
        else // vertex conflict
        {
            if (!cardinal1)
                cardinal1 = paths[get<0>(*con)]->at(get<4>(*con)).width==1;
            if (!cardinal2)
                cardinal2 = paths[get<1>(*con)]->at(get<4>(*con)).width==1;
        }
        if (cardinal1 && cardinal2)
        {
            parent.cardinalConf.push_back(con);
            if (h_type == heuristics_type::NONE)
                return;
            continue;
        }
        else if (rectangleReasoning // rectangle reasoning using MDDs
            && get<3>(*con) < 0 // vertex conflict
            && (int)paths[get<0>(*con)]->size() > get<4>(*con) && (int)paths[get<1>(*con)]->size() > get<4>(*con)) // conflict happens before agents reach their goals
        {
            //Rectangle reasoning for semi and non cardinal vertex conflicts
            int a1 = get<0>(*con);
            int a2 = get<1>(*con);
            int timestep = get<4>(*con);
            std::list<int>    s1s = getStartCandidates(*paths[a1], timestep, ml->cols);
            std::list<int>    g1s = getGoalCandidates(*paths[a1], timestep, ml->cols);
            std::list<int>    s2s = getStartCandidates(*paths[a2], timestep, ml->cols);
            std::list<int>    g2s = getGoalCandidates(*paths[a2], timestep, ml->cols);

            // Try all possible combinations
            bool found = false;
            std::shared_ptr<tuple<int, int, int, int, int>> conflict;
            int type = -1;
            int area = 0;
            for (int t1_start : s1s)
            {
                for (int t1_end : g1s)
                {
                    int s1 = paths[a1]->at(t1_start).location;
                    int g1 = paths[a1]->at(t1_end).location;
                    if (!isManhattanOptimal(s1, g1, t1_end - t1_start, ml->cols))
                        continue;
                    for (int t2_start : s2s)
                    {
                        for (int t2_end : g2s)
                        {
                            int s2 = paths[a2]->at(t2_start).location;
                            int g2 = paths[a2]->at(t2_end).location;
                            if (!isManhattanOptimal(s2, g2, t2_end - t2_start, ml->cols))
                                continue;
                            if (!isRectangleConflict(s1, s2, g1, g2, ml->cols))
                                continue;
                            std::pair<int, int> Rg = getRg(std::make_pair(s1 / ml->cols, s1 %  ml->cols), std::make_pair(g1 / ml->cols, g1 %  ml->cols),
                                std::make_pair(g2 / ml->cols, g2 %  ml->cols));
                            std::pair<int, int> Rs = getRs(std::make_pair(s1 / ml->cols, s1 %  ml->cols), std::make_pair(s2 / ml->cols, s2 %  ml->cols),
                                std::make_pair(g1 / ml->cols, g1 %  ml->cols));
                            int new_area = (abs(Rs.first - Rg.first) + 1) * (abs(Rs.second - Rg.second) + 1);
                            int new_type = classifyRectangleConflict(s1, s2, g1, g2, Rg, ml->cols);
                            if (new_type > type || (new_type == type && new_area > area))
                            {
                                conflict = std::shared_ptr<tuple<int, int, int, int, int>>
                                    (new tuple<int, int, int, int, int>(get<0>(*con), get<1>(*con), -1 - Rg.first *  ml->cols - Rg.second, t1_start, t2_start));
                                type = new_type;
                                area = new_area;
                            }
                        }
                    }
                }
            }
            if (type == 2)
            {
                parent.cardinalConf.push_back(conflict);
                continue;
            }
            else if (type == 1 && !findRectangleConflict(parent.parent, *conflict))
            {
                parent.rectSemiConf.push_back(conflict);
            }
            else if (type == 0 && !findRectangleConflict(parent.parent, *conflict))
            {
                parent.rectNonConf.push_back(conflict);
            }
        }
        if (cardinal1 || cardinal2)
        {
            parent.semiConf.push_back(con);
        }
        else
        {
            parent.nonConf.push_back(con);
        }
    }

    // remove conflicts that cannot be chosen, to save some memory
    list<pair<int, int>> highPriorityPairs;
    /*
    removeLowPriorityConflicts(parent.cardinalConf, highPriorityPairs, paths, ml->cols);
    removeLowPriorityConflicts(parent.rectSemiConf, highPriorityPairs, paths, ml->cols);
    removeLowPriorityConflicts(parent.semiConf, highPriorityPairs, paths, ml->cols);
    removeLowPriorityConflicts(parent.rectNonConf, highPriorityPairs, paths, ml->cols);
    removeLowPriorityConflicts(parent.nonConf, highPriorityPairs, paths, ml->cols);
    */
    
    while(!parent.allConf.empty())parent.allConf.pop_front();
    //while(!parent.confType.empty())parent.confType.pop_front();
    parent.confType.clear();
    for(auto c:parent.cardinalConf)parent.allConf.push_back(c),parent.confType.push_back(1);
    for(auto c:parent.rectSemiConf)parent.allConf.push_back(c),parent.confType.push_back(0);
    for(auto c:parent.semiConf)parent.allConf.push_back(c),parent.confType.push_back(0);
    for(auto c:parent.rectNonConf)parent.allConf.push_back(c),parent.confType.push_back(-1);
    for(auto c:parent.nonConf)parent.allConf.push_back(c),parent.confType.push_back(-1);
}

void ICBSSearch::removeLowPriorityConflicts(std::list<std::shared_ptr<Conflict>>& conflicts,
    list<pair<int, int>>& highPriorityPairs, const vector<vector<PathEntry>*>& paths, int num_col)
{
    if(conflicts.empty())
        return;
    list<pair<int, int>> got;
    for (auto p1 = conflicts.begin(); p1 != conflicts.end();)
    {
        int a1 = get<0>(**p1), a2 = get<1>(**p1);
        pair<int, int> ag(min(a1, a2), max(a1, a2));
        bool found = false;
        for (auto p : highPriorityPairs)
        {
            if (p == ag)
            {
                found = true;
                p1 = conflicts.erase(p1);
                break;
            }
        }
        if(found)
            continue;

        auto p2 = p1;
        p2++;
        bool keepP1 = true;
        while (p2 != conflicts.end())
        {
            if ((get<0>(**p2) == a1 && get<1>(**p2) == a2) || (get<1>(**p2) == a1 && get<0>(**p2) == a2))
            {
                if (max(get<4>(**p1), get<4>(**p2)) >= (int)min(paths[a1]->size(), paths[a2]->size()))
                    keepP1 = get<4>(**p1) <= get<4>(**p2);
                else
                {
                    int time1, time2;
                    if (get<2>(**p1) < 0) // rectangle conflict
                        time1 = getRectangleTime(**p1, paths, num_col);
                    else // vertex/edge conflict
                        time1 = std::get<4>(**p1);
                    if (get<2>(**p2) < 0) // rectangle conflict
                        time2 = getRectangleTime(**p2, paths, num_col);
                    else // vertex/edge conflict
                        time2 = std::get<4>(**p2);
                    keepP1 = time1 <= time2;
                }
                if (keepP1) //delete p2
                    p2 = conflicts.erase(p2);
                else //delete p1
                {
                    p1 = conflicts.erase(p1);
                    break;
                }
            }
            else
                ++p2;
        }
        if (keepP1)
            ++p1;
    }

    for (auto conflict : conflicts)
    {
        int a1 = get<0>(*conflict), a2 = get<1>(*conflict);
        highPriorityPairs.push_back(make_pair(min(a1, a2), max(a1, a2)));
    }
}

bool ICBSSearch::findPathForSingleAgent(ICBSNode*  node, int ag, int lowerbound)
{
    // extract all constraints on agent ag
    ICBSNode* curr = node;
    vector < list< pair<int, int> > > cons_vec;
    int minLength = collectConstraints(curr, ag, cons_vec);
    
    // build reservation table
    CAT cat(node->makespan + 1);  // initialized to false
    updateReservationTable(cat, ag, *node);
    // find a path w.r.t cons_vec (and prioretize by res_table).
    bool foundSol = search_engines[ag]->findPath(node->path, cons_vec, cat, std::max(minLength, lowerbound));
    LL_num_expanded += search_engines[ag]->num_expanded;
    LL_num_generated += search_engines[ag]->num_generated;
    //delete cons_vec;
    if (foundSol)
    {
        node->g_val = node->g_val - paths[ag]->size() + node->path.size();
        paths[ag] = &node->path;
        node->makespan = std::max(node->makespan, node->path.size() - 1);
        return true;
    }
    else
    {
        return false;
    }
}

bool ICBSSearch::generateChild(ICBSNode*  node, ICBSNode* parent,bool choosing)
{
    node->parent = parent;
    node->g_val = parent->g_val;
    node->makespan = parent->makespan;
    node->depth = parent->depth + 1;

    std::clock_t t1;

    t1 = std::clock();


    int lowerbound;
    if (get<2>(*parent->conflict) < 0) // both constraints are at goal locations
        lowerbound = 0; // lowver bound here is meaningless
    else if (get<4>(*parent->conflict) >= (int)paths[node->agent_id]->size()) //conflict happens after agent reaches its goal
        lowerbound = get<4>(*parent->conflict) + 1;
    //else if (parent->h_val > 0) // the chosen conflict is cardinal
    //    lowerbound = (int)paths[node->agent_id]->size();
    else
        lowerbound = (int)paths[node->agent_id]->size() - 1;
        
    if (!findPathForSingleAgent(node, node->agent_id, lowerbound))
        return false;

    
    runtime_lowlevel += (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;
    
    //Estimate h value
    node->h_val = 0;
    if (parent->h_val  == 0);
    else if (parent->conflictGraph.empty())
    {
        node->h_val = parent->h_val - 1; // stronger pathmax
    }
    else
    {
        int maxWeight = 0;
        boost::unordered_map<int, int>::iterator got;
        for (auto e : parent->conflictGraph)
        {
            if ((e.first / num_of_agents == node->agent_id || e.first % num_of_agents == node->agent_id) && e.second > maxWeight)
            {
                maxWeight = e.second;
                if (maxWeight >= parent->h_val)
                    break;
            }
        }
        if (maxWeight < parent->h_val)
            node->h_val = parent->h_val - maxWeight; // stronger pathmax
    }
    node->h_val = std::max(node->h_val, parent->f_val - node->g_val); // pathmax
    node->f_val = node->g_val + node->h_val;

    t1 = std::clock();
    findConflicts(*node);
    runtime_conflictdetection += (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;
    copyConflictGraph(*node, *node->parent);

    if(!choosing)
    {
        // update handles
        node->open_handle = open_list.push(node);
        HL_num_generated++;
        node->time_generated = HL_num_generated;
        if (node->f_val <= focal_list_threshold)
            node->focal_handle = focal_list.push(node);
        allNodes_table.push_back(node);
    }
    

    return true;
}

void ICBSSearch::copyConflictGraph(ICBSNode& child, const ICBSNode& parent)
{
    //copy conflict graph
    if (h_type == heuristics_type::DG || h_type == heuristics_type::WDG)
    {
        int a1 = child.agent_id; // the replaned agent
        int a2 = (a1 == get<0>(*parent.conflict)) ? get<1>(*parent.conflict) : get<0>(*parent.conflict); // the other agent in the chosen conflict
        int deltaCost = child.g_val - parent.g_val; // increased cost of agent a1
        for (auto e : parent.conflictGraph)
        {
            if (e.first / num_of_agents != a1 && e.first % num_of_agents != a1)
                child.conflictGraph[e.first] = e.second;
        }
    }
}

void ICBSSearch::printPaths() const
{
    for (int i = 0; i < num_of_agents; i++)
    {
        std::cout << "Agent " << i << " (" << paths_found_initially[i].size() - 1 << " -->" <<
            paths[i]->size() - 1 << "): ";
        for (size_t t = 0; t < paths[i]->size(); t++)
            std::cout << "(" << paths[i]->at(t).location / search_engines[0]->num_col << "," <<
                paths[i]->at(t).location % search_engines[0]->num_col << ")->";
        std::cout << std::endl;
    }
}


// adding new nodes to FOCAL (those with min-f-val*f_weight between the old and new LB)
void ICBSSearch::updateFocalList()
{
    ICBSNode* open_head = open_list.top();
    if (open_head->f_val > min_f_val)
    {
        if (screen == 3)
        {
            cout << "  Note -- FOCAL UPDATE!! from |FOCAL|=" << focal_list.size() << " with |OPEN|=" << open_list.size() << " to |FOCAL|=";
        }
        min_f_val = open_head->f_val;
        double new_focal_list_threshold = min_f_val * focal_w;
        for (ICBSNode* n : open_list)
        {
            if (n->f_val > focal_list_threshold &&
                n->f_val <= new_focal_list_threshold)
                n->focal_handle = focal_list.push(n);
        }
        focal_list_threshold = new_focal_list_threshold;
        if (screen == 3)
        {
            cout << focal_list.size() << endl;
        }
    }
}

void ICBSSearch::updateReservationTable(CAT& cat, int exclude_agent, const ICBSNode &node)
{
    for (int ag = 0; ag < num_of_agents; ag++)
    {
        if (ag != exclude_agent && paths[ag] != NULL)
        {
            for (size_t timestep = 0; timestep < node.makespan + 1; timestep++)
            {
                if (timestep >= paths[ag]->size())
                {
                    cat[timestep].insert(paths[ag]->back().location);
                }
                else// otherwise, return its location for that timestep
                {
                    int id = paths[ag]->at(timestep).location;
                    cat[timestep].insert(id);
                    if (timestep > 0 && paths[ag]->at(timestep - 1).location != id)
                    {
                        int prev_id = paths[ag]->at(timestep - 1).location;
                        cat[timestep].insert((1 + id) * ml->cols * ml->rows + prev_id);
                    }
                }
            }
        }
    }
}

void ICBSSearch::printResults() const
{
    if(solution_cost >= 0) // solved
        cout << "Optimal,";
    else if(solution_cost == -1) // time_out
        cout << "Timeout,";
    else if(solution_cost == -2) // no solution
        cout << "No solutions,";
    else if (solution_cost == -3) // nodes out
        cout << "Nodesout,";

    std::cout << runtime << "," <<
        HL_num_expanded << "," << HL_num_generated << "," <<
        LL_num_expanded << "," << LL_num_generated << "," <<
        solution_cost << "," << min_f_val << "," <<
        dummy_start->g_val << "," << dummy_start->f_val << "," <<
        mvc_runtime << "," <<
        cardinal_pair_runtime << "," << cardinal_pair_num << "," <<
        cardinal_heuristic_num << "," << cardinal_heuristic_value << "," <<
        not_cardinal_pair_runtime << "," << not_cardinal_pair_num << "," <<
        not_cardinal_heuristic_num << "," << not_cardinal_heuristic_value << "," <<
        bookingHitTimes << "," << bookingSearchtime << "," << HL_num_heuristics << "," <<
        build_mdds_runtime << "," <<mvc_runtime<<","<< sameScore<<","<<diffScore<<","<<std::endl;
}

void ICBSSearch::saveResults(const std::string &fileName, const std::string &instanceName) const
{
    ofstream stats;
    stats.open(fileName, ios::app);
    stats << runtime << "," <<
        HL_num_expanded << "," << HL_num_generated << "," <<totalConflict/HL_num_expanded<<","<<
        LL_num_expanded << "," << LL_num_generated << "," <<
        solution_cost << "," << min_f_val << "," <<
        dummy_start->g_val << "," << dummy_start->f_val << "," <<
        mvc_runtime << "," <<
        cardinal_pair_runtime << "," <<
        cardinal_pair_num << "," <<
        cardinal_heuristic_num << "," << cardinal_heuristic_value << "," <<
        not_cardinal_pair_runtime << "," <<
        not_cardinal_pair_num << "," <<
        not_cardinal_heuristic_num << "," << not_cardinal_heuristic_value << "," <<
        bookingHitTimes << ", " << bookingSearchtime << "," << HL_num_heuristics << "," <<
        build_mdds_runtime << "," <<
        h_type << "," << PC << "," <<  chooseConflictTime<<","<<instanceName <<","<< sameScore<<","<<diffScore<<","<< endl;
    stats.close();
}

void ICBSSearch::printConflicts(const ICBSNode &curr) const
{
    for (list<std::shared_ptr<tuple<int, int, int, int, int>>>::const_iterator it = curr.cardinalConf.begin(); it != curr.cardinalConf.end(); ++it)
    {
        std::cout << "Cardinal " << **it << std::endl;
    }
    for (list<std::shared_ptr<tuple<int, int, int, int, int>>>::const_iterator it = curr.semiConf.begin(); it != curr.semiConf.end(); ++it)
    {
        std::cout << "Semi " << **it << std::endl;
    }
    for (list<std::shared_ptr<tuple<int, int, int, int, int>>>::const_iterator it = curr.nonConf.begin(); it != curr.nonConf.end(); ++it)
    {
        std::cout << "Non " << **it << std::endl;
    }
    for (list<std::shared_ptr<tuple<int, int, int, int, int>>>::const_iterator it = curr.unknownConf.begin(); it != curr.unknownConf.end(); ++it)
    {
        std::cout << "Unknown " << **it << std::endl;
    }
}


void ICBSSearch::printStrategy() const
{
    switch (h_type)
    {
    case heuristics_type::NONE:
        if(PC)
            std::cout << "    ICBS: ";
        else
            std::cout << "     CBS: ";
        break;
    case heuristics_type::CG:
        std::cout << " CBSH+CG: ";
        break;
    case heuristics_type::DG:
        std::cout << " CBSH+DG: ";
        break;
    case heuristics_type::WDG:
        std::cout << "CBSH+WDG: ";
        break;
    default:
        exit(10);
    }
}

bool ICBSSearch::runICBSSearch()
{
    if(screen > 0) // 1 or 2
        printStrategy();
    // set timer
    start = std::clock();
    std::clock_t t1;

    if(useoracle==4){
       svm_rank.resize(100,0);
       freopen(svm_file.c_str(),"r",stdin);
       int id;double weight;
       while(scanf("%d %lf",&id,&weight)!=EOF)svm_rank[id-1]=weight;
       fclose(stdin);
    }
       
    // start is already in the open_list
    while (!focal_list.empty() && !solution_found)
    {
        if (min_f_val >= cost_upperbound)
        {
            solution_cost = min_f_val;
            solution_found = false;
            break;
        }
        runtime = (std::clock() - start) * 1000 / CLOCKS_PER_SEC;
        if (runtime > time_limit)
        {  // timeout
            solution_cost = -1;
            solution_found = false;
            break;
        }
        ICBSNode* curr = focal_list.top();
        focal_list.pop();
        open_list.erase(curr->open_handle);
        // takes the paths_found_initially and UPDATE all constrained paths found for agents from curr to dummy_start (and lower-bounds)
        t1 = std::clock();
        updatePaths(curr);
        runtime_updatepaths += (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;
        
        if (curr->num_of_collisions == 0) //no conflicts
        {// found a solution (and finish the while look)
            solution_found = true;
            solution_cost = curr->g_val;
            goal_node = curr;
            break;
        }
        else if (h_type == heuristics_type::NONE) // No heuristics
        {
            t1 = std::clock();
            if(PC) // priortize conflicts
                classifyConflicts(*curr);
            cerr<<"wrong"<<endl;
            //chooseConflict2(*curr);
            //chooseConflict2(*curr);
            runtime_conflictdetection += (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;
        }
        else if(curr->conflict == NULL) //use h value, and h value has not been computed yet
        {
            if (screen == 3)
            {
                std::cout << std::endl << "****** Compute h for #" << curr->time_generated << " with f= " << curr->g_val <<
                    "+" << curr->h_val << " (";
                for (int i = 0; i < num_of_agents; i++)
                    std::cout << paths[i]->size() - 1 << ", ";
                std::cout << ") and #conflicts = " << curr->num_of_collisions << std::endl;
            }

            t1 = std::clock();
            if (PC) // priortize conflicts
                classifyConflicts(*curr);
            runtime_conflictdetection += (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;

            t1 = std::clock();
            int h = computeHeuristics(*curr);
            double runtime_h = (std::clock() - t1) * 1000.0 / CLOCKS_PER_SEC;
            runtime_computeh += runtime_h;
            HL_num_heuristics++;

            if (h < 0) // no solution, so prune this node
            {
                curr->clear();
                if (open_list.size() == 0)
                {
                    solution_found = false;
                    break;
                }
                updateFocalList();
                continue;
            }

            curr->h_val = std::max(h, curr->h_val); // use consistent h values
            curr->f_val = curr->g_val + curr->h_val;

            if(screen == 2)
                curr->printConflictGraph(num_of_agents);

            //if (useoracle==0 or rand()%10==0)chooseConflict(*curr);
            if (useoracle==0||(useoracle==4&&curr->cardinalConf.size()>1))chooseConflict(*curr);
            //if (useoracle==0||(useoracle==4&&HL_num_expanded%2==0&&curr->cardinalConf.size()>=1))chooseConflict(*curr); // model 4
            //if (useoracle==0||(useoracle==4&&HL_num_expanded%2==0&&curr->cardinalConf.size()>1))chooseConflict(*curr); //model 5
            else if (useoracle==4)chooseConflict3(*curr);//large
            //else if (useoracle==4)chooseConflict3_Large(*curr);//large
            //else if (useoracle==4)chooseConflict3_Old(*curr);//feature5
            else  chooseConflict2(*curr);
            if (curr->f_val > focal_list_threshold)
            {
                if (screen == 3)
                {
                    std::cout << "Reinsert the node with f =" << curr->g_val << "+" << curr->h_val << std::endl;
                }

                curr->open_handle = open_list.push(curr);
                updateFocalList();
                continue;
            }
        }


         //Expand the node
        HL_num_expanded++;
        if(HL_num_expanded%1000==0)cerr<<HL_num_expanded<<" ";
        totalConflict+=curr->num_of_collisions;
        curr->time_expanded = HL_num_expanded;
        if(screen == 3)
            std::cout << "Expand Node " << curr->time_generated << " ( " << curr->f_val << "= " << curr->g_val << " + " <<
                curr->h_val << " ) on conflict " << *curr->conflict << std::endl;
        ICBSNode* n1 = new ICBSNode();
        ICBSNode* n2 = new ICBSNode();
        
        n1->agent_id = get<0>(*curr->conflict);
        n2->agent_id = get<1>(*curr->conflict);

        if (get<2>(*curr->conflict) < 0) // Rectangle conflict
        {
            int Rg = -1 - get<2>(*curr->conflict);
            int S1_t = get<3>(*curr->conflict);
            int S2_t = get<4>(*curr->conflict);
            const MDD* mdd1 = buildMDD(*curr, n1->agent_id);
            const MDD* mdd2 = buildMDD(*curr, n2->agent_id);
            addModifiedBarrierConstraints(*paths[get<0>(*curr->conflict)], *paths[get<1>(*curr->conflict)],
                mdd1, mdd2, S1_t, S2_t, Rg, ml->cols, n1->constraints, n2->constraints);
        }
        else if (get<3>(*curr->conflict) < 0) // vertex conflict
        {
            n1->constraints.push_back(make_tuple(get<2>(*curr->conflict), -1, get<4>(*curr->conflict)));
            n2->constraints.push_back(make_tuple(get<2>(*curr->conflict), -1, get<4>(*curr->conflict)));
            countResolveforNode[get<2>(*curr->conflict)]++;
        }
        else // edge conflict
        {
            n1->constraints.push_back(make_tuple(get<2>(*curr->conflict), get<3>(*curr->conflict), get<4>(*curr->conflict)));
            n2->constraints.push_back(make_tuple(get<3>(*curr->conflict), get<2>(*curr->conflict), get<4>(*curr->conflict)));
            countResolveforNode[get<2>(*curr->conflict)]+=0.5;
            countResolveforNode[get<3>(*curr->conflict)]+=0.5;
        }

        countResolveforAgent[n1->agent_id]++;
        countResolveforAgent[n2->agent_id]++;
        bool Sol1 = false, Sol2 = false;
        vector<vector<PathEntry>*> copy(paths);
        Sol1 = generateChild(n1, curr, 0);

        if (screen == 3 && Sol1)
        {
            std::cout << "Generate #" << n1->time_generated
                << " with cost " << n1->g_val
                << " and " << n1->num_of_collisions << " conflicts " << std::endl;
        }
        if (Sol1 && n1->f_val == min_f_val && n1->num_of_collisions == 0) //no conflicts
        {// found a solution (and finish the while look)
            solution_found = true;
            solution_cost = n1->g_val;
            goal_node = n1;
            break;
        }
        else if(!Sol1)
        {
            countResolveforAgent[n1->agent_id]--;
            delete (n1);
            n1 = NULL;
        }
        paths = copy;
        

        Sol2 = generateChild(n2, curr, 0);
        if (screen == 3 && Sol2)
        {
            std::cout << "Generate #" << n2->time_generated
                << " with cost " << n2->g_val
                << " and " << n2->num_of_collisions << " conflicts " << std::endl;

        }
        if (Sol2 && n2->f_val == min_f_val && n2->num_of_collisions == 0) //no conflicts
        {// found a solution (and finish the while look)
            solution_found = true;
            solution_cost = n2->g_val;
            goal_node = n2;
            break;
        }
        else if(!Sol2)
        {
            countResolveforAgent[n2->agent_id]--;
            delete (n2);
            n2 = NULL;
        }

        curr->clear();

        releaseMDDMemory(curr->agent_id);

        if (open_list.size() == 0)
        {
            solution_found = false;
            break;
        }
        updateFocalList();

    }  // end of while loop


    runtime = (std::clock() - start) * 1000.0 / CLOCKS_PER_SEC;
    if (solution_found && !validateSolution())
    {
        std::cout << "Solution invalid!!!" << std::endl;
        printPaths();
        exit(-1);
    }
    //cerr<<screen<<endl;
    if (screen > 0) // 1 or 2
        printResults();
    return solution_found;
}


void ICBSSearch::chooseConflict3_Large(ICBSNode &node)
{
    int currentmakespan=0;
    int currentsumofcost=0;
    int currentdepth=node.depth;
    std::vector<int> conflictsInvolved(num_of_agents,0);
    int numofConflict=node.num_of_collisions;
    int curr_fval=node.f_val;
    int curr_hval=node.h_val;
    int curr_gval=node.g_val;
    currentmakespan=node.makespan;
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        int a1 = get<0>(*conflict);
        int a2 = get<1>(*conflict);
        conflictsInvolved[a1]++;
        conflictsInvolved[a2]++;
    }
    std::vector<std::vector<double> > featureSet;
    std::vector<double> featureScore;
    featureSet.clear();featureScore.clear();
    clock_t t1=std::clock();
    int maxScore=-1;
    int minScore=1e8;
    std::shared_ptr<Conflict> chosenConflict;
    int numC=-1;
    
    freopen("inRunFeature.txt","w",stdout);//open output file

    int featureSize=0;
    updatePaths(&node);
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        node.conflict=conflict;
        //updatePaths(&node);
        numC++;
        
        int score=0;
            
        int K=5;
        //featuresCollected++;
        
        std::vector<double> features;
        features.clear();
        
        //score
        //features.push_back(score);

        //current makespan (1)
        //features.push_back(currentmakespan);//1
        //current depth
        //features.push_back(currentdepth);//2
        //current f,g,h value (3)
        //features.push_back(curr_fval);//3
        //features.push_back(curr_gval);//4
        //features.push_back(curr_hval);//5
        
        //number of conflicts(1)
        //features.push_back(node.num_of_collisions);//6

        //conflict type (2)
        if (get<3>(*node.conflict) < 0) //vertex conflict
            features.push_back(1),features.push_back(0);//7,8
        else features.push_back(0),features.push_back(1);
        
        //cardinality of the conflict (3)

        int conflictType=node.confType[numC];

        int conflicta=0,conflictb=0,conflictc=0;
        if(conflictType==1)conflicta=1;
        else if (conflictType==0)conflictb=1;
        else conflictc=1;
        features.push_back(conflicta);//9
        features.push_back(conflictb);//10
        features.push_back(conflictc);//11
        
        double pathL0=paths[get<0>(*node.conflict)]->size();
        double pathL1=paths[get<1>(*node.conflict)]->size();
        
        int u=get<2>(*node.conflict);
        int v=get<3>(*node.conflict);
        if (v==-1)v=u;
        
        //number of conflicts i,j resolved so far (2)
        features.push_back(countResolveforAgent[get<0>(*node.conflict)]); //12 x
        features.push_back(countResolveforAgent[get<1>(*node.conflict)]); //13 x
        
        features.push_back(countResolveforAgent[get<0>(*node.conflict)]+countResolveforAgent[get<1>(*node.conflict)]);//2nd batch to add
        
        //number of conflicts node u,v resolved so far (2)
        features.push_back(countResolveforNode[u]); //14
        features.push_back(countResolveforNode[v]); //15
        
        features.push_back(countResolveforNode[u]+countResolveforNode[v]);//new to add


        //number of conflicts i,j involved (2)
        features.push_back(conflictsInvolved[get<0>(*node.conflict)]); //16 x
        features.push_back(conflictsInvolved[get<1>(*node.conflict)]); //17 x
        features.push_back(conflictsInvolved[get<0>(*node.conflict)]+conflictsInvolved[get<1>(*node.conflict)]);//2nd batch to add

        //time of conflict (1)
        int timestep=get<4>(*node.conflict);
        features.push_back(timestep);//18
        
        //time of conflict/ currentmakespan (1)
        features.push_back(timestep*1.0/currentmakespan);//19
        
        //pathlength of agents i,j (2)
        features.push_back(pathL0);//20 ??
        features.push_back(pathL1);//21 ??
        features.push_back(pathL0+pathL1);//2nd batch to add
        features.push_back(abs(pathL0-pathL1));//2nd batch to add
        features.push_back(max(pathL0*1./pathL1,pathL1*1./pathL0));//2nd batch to add
        
        features.push_back(pathL0-shortestPathLen[get<0>(*node.conflict)]);//2nd batch to add
        features.push_back(pathL1-shortestPathLen[get<1>(*node.conflict)]);//2nd batch to add
        
        features.push_back(pathL0*1./shortestPathLen[get<0>(*node.conflict)]);//2nd batch to add
        features.push_back(pathL1*1./shortestPathLen[get<1>(*node.conflict)]);//2nd batch to add

        
        //pathlength / curr_fval (2)
        features.push_back(pathL0/curr_fval);//22
        features.push_back(pathL1/curr_fval);//23
        
        //whether i,j have reached goal (2)
        features.push_back(int(paths[get<0>(*node.conflict)]->size()>timestep));//24
        features.push_back(int(paths[get<1>(*node.conflict)]->size()>timestep));//25
        
        //pathlength-time (2)
        features.push_back(pathL0-timestep);//26
        features.push_back(pathL1-timestep);//27
        
        //(pathlength-time)/pathlength (2)
        features.push_back((pathL0-timestep)/pathL0);//28
        features.push_back((pathL1-timestep)/pathL1);//29
        
        //cerr<<"ok1"<<endl;
        
        //
        std::vector<int> nodeKstep(K,0);// nodes within k steps in the graph
        
        std::vector<int> distu;
        ml->findKNeighbor(u,v,K,distu,nodeKstep);
        //std::vector<int> distv=v==-1?distu:findKNeighbor(v,K);
        
        //cerr<<"ok2"<<endl;

        
        std::vector<int> conflictKstep(K+1,0);//conflicts within k steps in time-expanded graph
        std::vector<int> conflictKstep2(K+1,0);//conflicts within k steps in the graph
        std::vector<int> agentKstep(K+1,0);//agents within k steps in time-expanded graph
        
        for (std::shared_ptr<Conflict> tpCon:node.allConf)
        {
            int tu=get<2>(*tpCon);
            int tv=get<3>(*tpCon);
            int tt=get<4>(*tpCon);
            if (tv==-1)tv=tu;
            int tdist=min(distu[tu],distu[tv]);
            //int tdistv=min(distv[tu],distv[tv]);
            //int tdist=min(tdistu,tdistv);
            if(tdist>=1e7)continue;
            if(tdist+abs(tt-timestep)<=K)
                conflictKstep[tdist+abs(tt-get<4>(*node.conflict))]+=1;
            if(tdist<=K)conflictKstep2[tdist]+=1;
        }
        
        //cerr<<"ok3"<<endl;

        for(int a=0;a<num_of_agents;a++)
        {
            //int loc_a=paths[a]->at(timestep).location;
            int loc_a;
            int len_a=paths[a]->size();
            //if (distu[loc_a]>1e7 && distv[loc_a]>1e7) continue;
            for(int tt=timestep;tt>=0&&tt>=timestep-K;tt--)
            {
                if (tt<len_a) loc_a=paths[a]->at(tt).location;
                else loc_a=paths[a]->at(len_a-1).location;
                    
                if(timestep-tt+distu[loc_a]<=K)
                    agentKstep[timestep-tt+distu[loc_a]]+=1;
                else break;
            }
            for(int tt=timestep+1;tt<=timestep+K;tt++)
            {
                if (tt<len_a) loc_a=paths[a]->at(tt).location;
                else loc_a=paths[a]->at(len_a-1).location;
                if(tt-timestep+distu[loc_a]<=K)
                    agentKstep[tt-timestep+distu[loc_a]]+=1;
                else break;
            }
        }
        //cerr<<"ok4"<<endl;

        
        // (4k+3)
        features.insert(features.end(), conflictKstep.begin(), conflictKstep.end()); //30-35
        features.insert(features.end(), conflictKstep2.begin(), conflictKstep2.end());//36-41
        features.insert(features.end(), agentKstep.begin(), agentKstep.end());//42-47
        features.insert(features.end(), nodeKstep.begin(), nodeKstep.end());//48-52
        featureSet.push_back(features);
        featureScore.push_back(score);
        featureSize=features.size();
        /*
        printf("0 qid:1 ");
        for(int i=1;i<=features.size();i++)
            printf("%d:%.4lf ",i,features[i-1]);
        puts("");
         */

            /*
            freopen(featureFile.c_str(),"a",stdout);
            for(auto feat:features)printf("%.4lf ",feat);
            puts("");
            fclose(stdout);*/
        
        //releaseMDDMemory(curr->agent_id);
    }
    //if (minScore==maxScore)sameScore+=1;
    //else
    //{
    //freopen("inRunFeature.txt","w",stdout);
    //sort(featureScore.begin(),featureScore.end());
    //reverse(featureScore.begin(),featureScore.end());
    /*int tp=featureScore.size()/10;
    int tpj=tp+1;
    for(;tpj<featureScore.size();tpj++)
        if(featureScore[tp]!=featureScore[tpj])break;
    //printf("%d %d ",tp,tpj);
    if(tpj*1.0/featureScore.size()>=0.2)
    {
        int tpk=tp;
        for(;tp>0;tp--)
            if(featureScore[tpk]!=featureScore[tp])break;
        //printf("%d\n",tp);
    }
    //printf("\n");
    diffScore+=1;
    featuresCollected++;
    for(auto feature:featureSet)
    {
        feature.insert(feature.begin(),featuresCollected);
        if(feature[1]+1e-8>=featureScore[tp])feature.insert(feature.begin(),1);
        else feature.insert(feature.begin(),0);
        for(auto feat:feature)
            if(fabs(feat-int(feat+0.01))<1e-8)printf("%d ",int(feat+0.01));
            else printf("%.4lf ",feat);
        puts("");
    }
    puts("");*/
    vector<double> norm(featureSize,0);
    for(auto features:featureSet)
    {
        for(int i=0;i<features.size();i++)
            norm[i]=max(norm[i],abs(features[i]));
    }
    for(int i=0;i<norm.size();i++)
        if(norm[i]<1e-8)norm[i]=1;
    for(auto features:featureSet)
    {
        printf("0 qid:1 ");
        for(int i=1;i<=features.size();i++)
            printf("%d:%.4lf ",i,features[i-1]/norm[i-1]);
        puts("");
    }
    fclose(stdout);
    //exit(0);
    //}
    
    // command ./svm_rank_classify inRunFeature.txt featureMax5/model inRunPrediction.txt
    system("./svm_rank_classify inRunFeature.txt featureMax5/model1/model inRunPrediction.txt");
    //system("./svm_rank_classify inRunFeature.txt featureLarge100/model inRunPrediction.txt");
    int currScore=-1e9;
    freopen("inRunPrediction.txt","r",stdin);
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        double score;
        scanf("%lf",&score);
        //cerr<<score<<" ";
        if(score>currScore)
        {
            currScore=score;
            chosenConflict=conflict;
        }
    }
    //cerr<<node.allConf.size()<<endl;
    
    fclose(stdin);
    //exit(0);
    node.conflict = chosenConflict;
    chooseConflictTime+=(std::clock() - t1) * 1000.0/ CLOCKS_PER_SEC;

}

void ICBSSearch::chooseConflict3_Old(ICBSNode &node)
{
    int currentmakespan=0;
    int currentsumofcost=0;
    int currentdepth=node.depth;
    std::vector<int> conflictsInvolved(num_of_agents,0);
    int numofConflict=node.num_of_collisions;
    int curr_fval=node.f_val;
    int curr_hval=node.h_val;
    int curr_gval=node.g_val;
    currentmakespan=node.makespan;
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        int a1 = get<0>(*conflict);
        int a2 = get<1>(*conflict);
        conflictsInvolved[a1]++;
        conflictsInvolved[a2]++;
    }
    std::vector<std::vector<double> > featureSet;
    std::vector<double> featureScore;
    featureSet.clear();featureScore.clear();
    clock_t t1=std::clock();
    int maxScore=-1;
    int minScore=1e8;
    std::shared_ptr<Conflict> chosenConflict;
    int numC=-1;
    
    freopen("inRunFeature.txt","w",stdout);//open output file

    int featureSize=0;
    updatePaths(&node);
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        node.conflict=conflict;
        //updatePaths(&node);
        numC++;
        
        int score=0;
            
        int K=5;
        //featuresCollected++;
        
        std::vector<double> features;
        features.clear();
        
        //score
        //features.push_back(score);

        //current makespan (1)
        features.push_back(currentmakespan);//1
        //current depth
        features.push_back(currentdepth);//2
        //current f,g,h value (3)
        features.push_back(curr_fval);//3
        features.push_back(curr_gval);//4
        features.push_back(curr_hval);//5
        
        //number of conflicts(1)
        features.push_back(node.num_of_collisions);//6

        //conflict type (2)
        if (get<3>(*node.conflict) < 0) //vertex conflict
            features.push_back(1),features.push_back(0);//7,8
        else features.push_back(0),features.push_back(1);
        
        //cardinality of the conflict (3)

        int conflictType=node.confType[numC];

        int conflicta=0,conflictb=0,conflictc=0;
        if(conflictType==1)conflicta=1;
        else if (conflictType==0)conflictb=1;
        else conflictc=1;
        features.push_back(conflicta);//9
        features.push_back(conflictb);//10
        features.push_back(conflictc);//11
        
        double pathL0=paths[get<0>(*node.conflict)]->size();
        double pathL1=paths[get<1>(*node.conflict)]->size();
        
        int u=get<2>(*node.conflict);
        int v=get<3>(*node.conflict);
        if (v==-1)v=u;
        
        //number of conflicts i,j resolved so far (2)
        features.push_back(countResolveforAgent[get<0>(*node.conflict)]); //12 x
        features.push_back(countResolveforAgent[get<1>(*node.conflict)]); //13 x
        
        //features.push_back(countResolveforAgent[get<0>(*node.conflict)]+countResolveforAgent[get<1>(*node.conflict)]);//new to add
        
        //number of conflicts node u,v resolved so far (2)
        features.push_back(countResolveforNode[u]); //14
        features.push_back(countResolveforNode[v]); //15
        
        //features.push_back(countResolveforNode[u]+countResolveforNode[v]);//new to add


        //number of conflicts i,j involved (2)
        features.push_back(conflictsInvolved[get<0>(*node.conflict)]); //16 x
        features.push_back(conflictsInvolved[get<1>(*node.conflict)]); //17 x
        //features.push_back(conflictsInvolved[get<0>(*node.conflict)]+conflictsInvolved[get<1>(*node.conflict)]);//new to add

        //time of conflict (1)
        int timestep=get<4>(*node.conflict);
        features.push_back(timestep);//18
        
        //time of conflict/ currentmakespan (1)
        features.push_back(timestep*1.0/currentmakespan);//19
        
        //pathlength of agents i,j (2)
        features.push_back(pathL0);//20 ??
        features.push_back(pathL1);//21 ??
        //features.push_back(pathL0+pathL1);
        //features.push_back(abs(pathL0-pathL1));
        //features.push_back(max(pathL0*1./pathL1,pathL1*1./pathL0));
        
        //features.push_back(pathL0-shortestPathLen[get<0>(*node.conflict)]);
        //features.push_back(pathL1-shortestPathLen[get<1>(*node.conflict)]);
        
        //features.push_back(pathL0*1./shortestPathLen[get<0>(*node.conflict)]);
        //features.push_back(pathL1*1./shortestPathLen[get<1>(*node.conflict)]);

        
        //pathlength / curr_fval (2)
        features.push_back(pathL0/curr_fval);//22
        features.push_back(pathL1/curr_fval);//23
        
        //whether i,j have reached goal (2)
        features.push_back(int(paths[get<0>(*node.conflict)]->size()>timestep));//24
        features.push_back(int(paths[get<1>(*node.conflict)]->size()>timestep));//25
        
        //pathlength-time (2)
        features.push_back(pathL0-timestep);//26
        features.push_back(pathL1-timestep);//27
        
        //(pathlength-time)/pathlength (2)
        features.push_back((pathL0-timestep)/pathL0);//28
        features.push_back((pathL1-timestep)/pathL1);//29
        
        //cerr<<"ok1"<<endl;
        
        //
        std::vector<int> nodeKstep(K,0);// nodes within k steps in the graph
        
        std::vector<int> distu;
        ml->findKNeighbor(u,v,K,distu,nodeKstep);
        //std::vector<int> distv=v==-1?distu:findKNeighbor(v,K);
        
        //cerr<<"ok2"<<endl;

        
        std::vector<int> conflictKstep(K+1,0);//conflicts within k steps in time-expanded graph
        std::vector<int> conflictKstep2(K+1,0);//conflicts within k steps in the graph
        std::vector<int> agentKstep(K+1,0);//agents within k steps in time-expanded graph
        
        for (std::shared_ptr<Conflict> tpCon:node.allConf)
        {
            int tu=get<2>(*tpCon);
            int tv=get<3>(*tpCon);
            int tt=get<4>(*tpCon);
            if (tv==-1)tv=tu;
            int tdist=min(distu[tu],distu[tv]);
            //int tdistv=min(distv[tu],distv[tv]);
            //int tdist=min(tdistu,tdistv);
            if(tdist>=1e7)continue;
            if(tdist+abs(tt-timestep)<=K)
                conflictKstep[tdist+abs(tt-get<4>(*node.conflict))]+=1;
            if(tdist<=K)conflictKstep2[tdist]+=1;
        }
        
        //cerr<<"ok3"<<endl;

        for(int a=0;a<num_of_agents;a++)
        {
            //int loc_a=paths[a]->at(timestep).location;
            int loc_a;
            int len_a=paths[a]->size();
            //if (distu[loc_a]>1e7 && distv[loc_a]>1e7) continue;
            for(int tt=timestep;tt>=0&&tt>=timestep-K;tt--)
            {
                if (tt<len_a) loc_a=paths[a]->at(tt).location;
                else loc_a=paths[a]->at(len_a-1).location;
                    
                if(timestep-tt+distu[loc_a]<=K)
                    agentKstep[timestep-tt+distu[loc_a]]+=1;
                else break;
            }
            for(int tt=timestep+1;tt<=timestep+K;tt++)
            {
                if (tt<len_a) loc_a=paths[a]->at(tt).location;
                else loc_a=paths[a]->at(len_a-1).location;
                if(tt-timestep+distu[loc_a]<=K)
                    agentKstep[tt-timestep+distu[loc_a]]+=1;
                else break;
            }
        }
        //cerr<<"ok4"<<endl;

        
        // (4k+3)
        features.insert(features.end(), conflictKstep.begin(), conflictKstep.end()); //30-35
        features.insert(features.end(), conflictKstep2.begin(), conflictKstep2.end());//36-41
        features.insert(features.end(), agentKstep.begin(), agentKstep.end());//42-47
        features.insert(features.end(), nodeKstep.begin(), nodeKstep.end());//48-52
        featureSet.push_back(features);
        featureScore.push_back(score);
        featureSize=features.size();
        /*
        printf("0 qid:1 ");
        for(int i=1;i<=features.size();i++)
            printf("%d:%.4lf ",i,features[i-1]);
        puts("");
         */

            /*
            freopen(featureFile.c_str(),"a",stdout);
            for(auto feat:features)printf("%.4lf ",feat);
            puts("");
            fclose(stdout);*/
        
        //releaseMDDMemory(curr->agent_id);
    }
    //if (minScore==maxScore)sameScore+=1;
    //else
    //{
    //freopen("inRunFeature.txt","w",stdout);
    //sort(featureScore.begin(),featureScore.end());
    //reverse(featureScore.begin(),featureScore.end());
    /*int tp=featureScore.size()/10;
    int tpj=tp+1;
    for(;tpj<featureScore.size();tpj++)
        if(featureScore[tp]!=featureScore[tpj])break;
    //printf("%d %d ",tp,tpj);
    if(tpj*1.0/featureScore.size()>=0.2)
    {
        int tpk=tp;
        for(;tp>0;tp--)
            if(featureScore[tpk]!=featureScore[tp])break;
        //printf("%d\n",tp);
    }
    //printf("\n");
    diffScore+=1;
    featuresCollected++;
    for(auto feature:featureSet)
    {
        feature.insert(feature.begin(),featuresCollected);
        if(feature[1]+1e-8>=featureScore[tp])feature.insert(feature.begin(),1);
        else feature.insert(feature.begin(),0);
        for(auto feat:feature)
            if(fabs(feat-int(feat+0.01))<1e-8)printf("%d ",int(feat+0.01));
            else printf("%.4lf ",feat);
        puts("");
    }
    puts("");*/
    vector<double> norm(featureSize,0);
    for(auto features:featureSet)
    {
        for(int i=0;i<features.size();i++)
            norm[i]=max(norm[i],abs(features[i]));
    }
    for(int i=0;i<norm.size();i++)
        //if(norm[i]<1e-8)norm[i]=1;
        norm[i]=1;
    for(auto features:featureSet)
    {
        printf("0 qid:1 ");
        for(int i=1;i<=features.size();i++)
            printf("%d:%.4lf ",i,features[i-1]/norm[i-1]);
        puts("");
    }
    fclose(stdout);
    //exit(0);
    //}
    
    // command ./svm_rank_classify inRunFeature.txt featureMax5/model inRunPrediction.txt
    system("./svm_rank_classify inRunFeature.txt featureMax5/model1/model inRunPrediction.txt");
    //system("./svm_rank_classify inRunFeature.txt featureLarge100/model inRunPrediction.txt");
    int currScore=-1e9;
    freopen("inRunPrediction.txt","r",stdin);
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        double score;
        scanf("%lf",&score);
        //cerr<<score<<" ";
        if(score>currScore)
        {
            currScore=score;
            chosenConflict=conflict;
        }
    }
    //cerr<<node.allConf.size()<<endl;
    
    fclose(stdin);
    //exit(0);
    node.conflict = chosenConflict;
    chooseConflictTime+=(std::clock() - t1) * 1000.0/ CLOCKS_PER_SEC;

}

void ICBSSearch::chooseConflict3(ICBSNode &node)
{
    int currentmakespan=0;
    int currentsumofcost=0;
    int currentdepth=node.depth;
    std::vector<int> conflictsInvolved(num_of_agents,0);
    int numofConflict=node.num_of_collisions;
    int curr_fval=node.f_val;
    int curr_hval=node.h_val;
    int curr_gval=node.g_val;
    currentmakespan=node.makespan;
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        int a1 = get<0>(*conflict);
        int a2 = get<1>(*conflict);
        conflictsInvolved[a1]++;
        conflictsInvolved[a2]++;
    }
    std::vector<std::vector<double> > featureSet;
    std::vector<double> featureScore;
    featureSet.clear();featureScore.clear();
    clock_t t1=std::clock();
    int maxScore=-1;
    int minScore=1e8;
    std::shared_ptr<Conflict> chosenConflict;
    int numC=-1;
    int numCard=node.cardinalConf.size();
    //if(numCard<=1)numCard=1e8;
    //freopen("inRunFeature.txt","w",stdout);//open output file

    int featureSize=0;
    updatePaths(&node);
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        node.conflict=conflict;
        //updatePaths(&node);
        numC++;
        //if(numCard>1&&node.confType[numC]!=1)break;
        
     
        int score=0;
            
        int K=5;
        //featuresCollected++;
        
        std::vector<double> features;
        features.clear();
        
        //score
        //features.push_back(score);

        //current makespan (1)
        //features.push_back(currentmakespan);//1
        //current depth
        //features.push_back(currentdepth);//2
        //current f,g,h value (3)
        //features.push_back(curr_fval);//3
        //features.push_back(curr_gval);//4
        //features.push_back(curr_hval);//5
        
        //number of conflicts(1)
        //features.push_back(node.num_of_collisions);//6

        //conflict type (2)
        if (get<3>(*node.conflict) < 0) //vertex conflict
            features.push_back(1),features.push_back(0);//7,8
        else features.push_back(0),features.push_back(1);
        
        //cardinality of the conflict (3)

        int conflictType=node.confType[numC];

        int conflicta=0,conflictb=0,conflictc=0;
        if(conflictType==1)conflicta=1;
        else if (conflictType==0)conflictb=1;
        else conflictc=1;
        features.push_back(conflicta);//9
        features.push_back(conflictb);//10
        features.push_back(conflictc);//11
        
        double pathL0=paths[get<0>(*node.conflict)]->size();
        double pathL1=paths[get<1>(*node.conflict)]->size();
        
        int u=get<2>(*node.conflict);
        int v=get<3>(*node.conflict);
        if (v==-1)v=u;
        
        //number of conflicts i,j resolved so far (2)
        features.push_back(countResolveforAgent[get<0>(*node.conflict)]); //12 x
        features.push_back(countResolveforAgent[get<1>(*node.conflict)]); //13 x
        
        features.push_back(countResolveforAgent[get<0>(*node.conflict)]+countResolveforAgent[get<1>(*node.conflict)]);//2nd batch to add
        
        //number of conflicts node u,v resolved so far (2)
        features.push_back(countResolveforNode[u]); //14
        features.push_back(countResolveforNode[v]); //15
        
        features.push_back(countResolveforNode[u]+countResolveforNode[v]);//new to add


        //number of conflicts i,j involved (2)
        features.push_back(conflictsInvolved[get<0>(*node.conflict)]); //16 x
        features.push_back(conflictsInvolved[get<1>(*node.conflict)]); //17 x
        features.push_back(conflictsInvolved[get<0>(*node.conflict)]+conflictsInvolved[get<1>(*node.conflict)]);//2nd batch to add

        //time of conflict (1)
        int timestep=get<4>(*node.conflict);
        features.push_back(timestep);//18
        
        //time of conflict/ currentmakespan (1)
        features.push_back(timestep*1.0/currentmakespan);//19
        
        //pathlength of agents i,j (2)
        features.push_back(pathL0);//20 ??
        features.push_back(pathL1);//21 ??
        features.push_back(pathL0+pathL1);//2nd batch to add
        features.push_back(abs(pathL0-pathL1));//2nd batch to add
        features.push_back(max(pathL0*1./pathL1,pathL1*1./pathL0));//2nd batch to add
        
        features.push_back(pathL0-shortestPathLen[get<0>(*node.conflict)]);//2nd batch to add
        features.push_back(pathL1-shortestPathLen[get<1>(*node.conflict)]);//2nd batch to add
        
        features.push_back(pathL0*1./shortestPathLen[get<0>(*node.conflict)]);//2nd batch to add
        features.push_back(pathL1*1./shortestPathLen[get<1>(*node.conflict)]);//2nd batch to add

        
        //pathlength / curr_fval (2)
        features.push_back(pathL0/curr_fval);//22
        features.push_back(pathL1/curr_fval);//23
        
        //whether i,j have reached goal (2)
        features.push_back(int(paths[get<0>(*node.conflict)]->size()>timestep));//24
        features.push_back(int(paths[get<1>(*node.conflict)]->size()>timestep));//25
        
        //pathlength-time (2)
        features.push_back(pathL0-timestep);//26
        features.push_back(pathL1-timestep);//27
        
        //(pathlength-time)/pathlength (2)
        features.push_back((pathL0-timestep)/pathL0);//28
        features.push_back((pathL1-timestep)/pathL1);//29
        
        //cerr<<"ok1"<<endl;
        
        //
        std::vector<int> nodeKstep(K,0);// nodes within k steps in the graph
        
        std::vector<int> distu;
        ml->findKNeighbor(u,v,K,distu,nodeKstep);
        //std::vector<int> distv=v==-1?distu:findKNeighbor(v,K);
        
        //cerr<<"ok2"<<endl;

        
        std::vector<int> conflictKstep(K+1,0);//conflicts within k steps in time-expanded graph
        std::vector<int> conflictKstep2(K+1,0);//conflicts within k steps in the graph
        std::vector<int> agentKstep(K+1,0);//agents within k steps in time-expanded graph
        
        for (std::shared_ptr<Conflict> tpCon:node.allConf)
        {
            int tu=get<2>(*tpCon);
            int tv=get<3>(*tpCon);
            int tt=get<4>(*tpCon);
            if (tv==-1)tv=tu;
            int tdist=min(distu[tu],distu[tv]);
            //int tdistv=min(distv[tu],distv[tv]);
            //int tdist=min(tdistu,tdistv);
            if(tdist>=1e7)continue;
            if(tdist+abs(tt-timestep)<=K)
                conflictKstep[tdist+abs(tt-get<4>(*node.conflict))]+=1;
            if(tdist<=K)conflictKstep2[tdist]+=1;
        }
        
        //cerr<<"ok3"<<endl;

        for(int a=0;a<num_of_agents;a++)
        {
            //int loc_a=paths[a]->at(timestep).location;
            int loc_a;
            int len_a=paths[a]->size();
            //if (distu[loc_a]>1e7 && distv[loc_a]>1e7) continue;
            for(int tt=timestep;tt>=0&&tt>=timestep-K;tt--)
            {
                if (tt<len_a) loc_a=paths[a]->at(tt).location;
                else loc_a=paths[a]->at(len_a-1).location;
                    
                if(timestep-tt+distu[loc_a]<=K)
                    agentKstep[timestep-tt+distu[loc_a]]+=1;
                else break;
            }
            for(int tt=timestep+1;tt<=timestep+K;tt++)
            {
                if (tt<len_a) loc_a=paths[a]->at(tt).location;
                else loc_a=paths[a]->at(len_a-1).location;
                if(tt-timestep+distu[loc_a]<=K)
                    agentKstep[tt-timestep+distu[loc_a]]+=1;
                else break;
            }
        }
        //cerr<<"ok4"<<endl;

        
        // (4k+3)
        features.insert(features.end(), conflictKstep.begin(), conflictKstep.end()); //30-35
        features.insert(features.end(), conflictKstep2.begin(), conflictKstep2.end());//36-41
        features.insert(features.end(), agentKstep.begin(), agentKstep.end());//42-47
        features.insert(features.end(), nodeKstep.begin(), nodeKstep.end());//48-52
        
        int a1 = min(get<0>(*node.conflict), get<1>(*node.conflict));
        int a2 = max(get<0>(*node.conflict), get<1>(*node.conflict));
        for(int tt=timestep-2;tt<=timestep+2;tt++)
        {
            int mddw1=0,mddw2=0;
            if(tt>=0)
            {
                if(tt<paths[a1]->size())mddw1=paths[a1]->at(tt).width;
                if(tt<paths[a2]->size())mddw2=paths[a2]->at(tt).width;
            }
            features.push_back(mddw1);
            features.push_back(mddw2);
        }//3rd new add
        
        int idx = a1 * num_of_agents + a2;
        auto got = node.conflictGraph.find(idx);
        if (got != node.conflictGraph.end())
            features.push_back(got->second);
        else {features.push_back(0);cerr<<"no weight in DG!"<<endl;exit(0);}
        
        
        featureSet.push_back(features);
        featureScore.push_back(score);
        featureSize=features.size();
        /*
        printf("0 qid:1 ");
        for(int i=1;i<=features.size();i++)
            printf("%d:%.4lf ",i,features[i-1]);
        puts("");
         */

            /*
            freopen(featureFile.c_str(),"a",stdout);
            for(auto feat:features)printf("%.4lf ",feat);
            puts("");
            fclose(stdout);*/
        
        //releaseMDDMemory(curr->agent_id);
    }
    //if (minScore==maxScore)sameScore+=1;
    //else
    //{
    //freopen("inRunFeature.txt","w",stdout);
    //sort(featureScore.begin(),featureScore.end());
    //reverse(featureScore.begin(),featureScore.end());
    /*int tp=featureScore.size()/10;
    int tpj=tp+1;
    for(;tpj<featureScore.size();tpj++)
        if(featureScore[tp]!=featureScore[tpj])break;
    //printf("%d %d ",tp,tpj);
    if(tpj*1.0/featureScore.size()>=0.2)
    {
        int tpk=tp;
        for(;tp>0;tp--)
            if(featureScore[tpk]!=featureScore[tp])break;
        //printf("%d\n",tp);
    }
    //printf("\n");
    diffScore+=1;
    featuresCollected++;
    for(auto feature:featureSet)
    {
        feature.insert(feature.begin(),featuresCollected);
        if(feature[1]+1e-8>=featureScore[tp])feature.insert(feature.begin(),1);
        else feature.insert(feature.begin(),0);
        for(auto feat:feature)
            if(fabs(feat-int(feat+0.01))<1e-8)printf("%d ",int(feat+0.01));
            else printf("%.4lf ",feat);
        puts("");
    }
    puts("");*/
    vector<double> norm(featureSize,0);
    vector<double> score_feat;
    score_feat.clear();
    for(auto features:featureSet)
    {
        for(int i=0;i<features.size();i++)
            norm[i]=max(norm[i],abs(features[i]));
    }
    for(int i=0;i<norm.size();i++)
        if(norm[i]<1e-8)norm[i]=1;
    for(auto features:featureSet)
    {
        //printf("0 qid:1 ");
        double score=0;
        for(int i=1;i<=features.size();i++)
            score+=svm_rank[i-1]*features[i-1]/norm[i-1];
            //printf("%d:%.4lf ",i,features[i-1]/norm[i-1]);
        //puts("");
        score_feat.push_back(score);
    }
    //fclose(stdout);
    //exit(0);
    //}
    
    // command ./svm_rank_classify inRunFeature.txt featureMax5/model inRunPrediction.txt
    //system("./svm_rank_classify inRunFeature.txt featureMax5/modelcombined inRunPrediction.txt");
    //system("./svm_rank_classify inRunFeature.txt featureMax/mode2 inRunPrediction.txt");//map 3 and map 5
    //system("./svm_rank_classify inRunFeature.txt featureLarge100/model inRunPrediction.txt");
    //system("./svm_rank_classify inRunFeature.txt featureLarge100/large_comb inRunPrediction.txt");
    
    
    int currScore=-1e9;numC=-1;
    //freopen("inRunPrediction.txt","r",stdin);
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        numC++;
        double score=score_feat[numC];
        
        //if(numC>=numCard)break;
        //if(numCard>1&&node.confType[numC]!=1)break;
        //if(scanf("%lf",&score)==EOF){cerr<<"EOF"<<endl;exit(0);};
        //cerr<<score<<" ";
        if(score>currScore)
        {
            currScore=score;
            chosenConflict=conflict;
        }
    }
    //cerr<<node.allConf.size()<<endl;
    
    //fclose(stdin);
    //exit(0);
    node.conflict = chosenConflict;
    chooseConflictTime+=(std::clock() - t1) * 1000.0/ CLOCKS_PER_SEC;

}


void ICBSSearch::chooseConflict2(ICBSNode &node)
{
    //cerr<<"wrong"<<endl;
    int currentmakespan=0;
    int currentsumofcost=0;
    int currentdepth=node.depth;
    std::vector<int> conflictsInvolved(num_of_agents,0);
    int numofConflict=node.num_of_collisions;
    int curr_fval=node.f_val;
    int curr_hval=node.h_val;
    int curr_gval=node.g_val;
    if (collectingFeatures)
    {
        currentmakespan=node.makespan;
        for (std::shared_ptr<Conflict> conflict: node.allConf)
        {
            int a1 = get<0>(*conflict);
            int a2 = get<1>(*conflict);
            conflictsInvolved[a1]++;
            conflictsInvolved[a2]++;
        }
    }
    std::vector<std::vector<double> > featureSet;
    std::vector<double> featureScore;
    featureSet.clear();featureScore.clear();
    clock_t t1=std::clock();
    double maxScore=-1;
    double minScore=1e8;
    
    std::shared_ptr<Conflict> chosenConflict;
    int numC=-1;
    for (std::shared_ptr<Conflict> conflict: node.allConf)
    {
        node.conflict=conflict;
        updatePaths(&node);
        numC++;
        ICBSNode* n1 = new ICBSNode();
        ICBSNode* n2 = new ICBSNode();
        
        n1->agent_id = get<0>(*node.conflict);
        n2->agent_id = get<1>(*node.conflict);
        
        int minConf=1e8;
        
        if (get<2>(*node.conflict) < 0) // Rectangle conflict
        {
            int Rg = -1 - get<2>(*node.conflict);
            int S1_t = get<3>(*node.conflict);
            int S2_t = get<4>(*node.conflict);
            const MDD* mdd1 = buildMDD(node, n1->agent_id);
            const MDD* mdd2 = buildMDD(node, n2->agent_id);
            addModifiedBarrierConstraints(*paths[get<0>(*node.conflict)], *paths[get<1>(*node.conflict)],
                mdd1, mdd2, S1_t, S2_t, Rg, ml->cols, n1->constraints, n2->constraints);
        }
        else if (get<3>(*node.conflict) < 0) // vertex conflict
        {
            n1->constraints.push_back(make_tuple(get<2>(*node.conflict), -1, get<4>(*node.conflict)));
            n2->constraints.push_back(make_tuple(get<2>(*node.conflict), -1, get<4>(*node.conflict)));
        }
        else // edge conflict
        {
            n1->constraints.push_back(make_tuple(get<2>(*node.conflict), get<3>(*node.conflict), get<4>(*node.conflict)));
            n2->constraints.push_back(make_tuple(get<3>(*node.conflict), get<2>(*node.conflict), get<4>(*node.conflict)));
        }

        bool Sol1 = false, Sol2 = false;
        vector<vector<PathEntry>*> copy(paths);
        Sol1 = generateChild(n1, &node,1);
        if (screen == 3 && Sol1)
        {
            std::cout << "Generate #" << n1->time_generated
                << " with cost " << n1->g_val
                << " and " << n1->num_of_collisions << " conflicts " << std::endl;
        }
        int h1=1e8,h2=1e8;
        
        int maxscore_1,maxscore_2;
        int minconf_1,minconf_2;
        
        if (Sol1 && n1->f_val == min_f_val && n1->num_of_collisions == 0) //no conflicts
        {// found a solution (and finish the while look)
            //solution_found = true;
            //solution_cost = n1->g_val;
            //goal_node = n1;
            //break;
            h1=n1->g_val;
        }
        else if(!Sol1)
        {
            delete (n1);
            n1 = NULL;
        }else if (useoracle==1 || useoracle==3)
        {
            updatePaths(n1);
            classifyConflicts(*n1);
            h1=computeHeuristics(*n1);
            n1->h_val = std::max(h1, n1->h_val); // use consistent h values
            n1->f_val = n1->g_val + n1->h_val;
            h1=n1->f_val;
            n1->clear();
        }//else h1=n1->num_of_collisions;
        if(n1!=NULL)minConf=min(minConf,n1->num_of_collisions);
        
        
        paths = copy;
        Sol2 = generateChild(n2, &node,1);
        if (screen == 3 && Sol2)
        {
            std::cout << "Generate #" << n2->time_generated
                << " with cost " << n2->g_val
                << " and " << n2->num_of_collisions << " conflicts " << std::endl;

        }
        if (Sol2 && n2->f_val == min_f_val && n2->num_of_collisions == 0) //no conflicts
        {// found a solution (and finish the while look)
            //solution_found = true;
            //solution_cost = n2->g_val;
            //goal_node = n2;
            h2=n2->g_val;
            //break;
        }
        else if(!Sol2)
        {
            delete (n2);
            n2 = NULL;
        }else if (useoracle==1 || useoracle==3)
        {
            updatePaths(n2);
            classifyConflicts(*n1);
            h2=computeHeuristics(*n2);
            n2->h_val = std::max(h2, n2->h_val); // use consistent h values
            n2->f_val = n2->g_val + n2->h_val;
            h2=n2->f_val;
            n2->clear();
        }//else h2=n2->num_of_collisions;
        if(n2!=NULL)minConf=min(minConf,n2->num_of_collisions);
        
        paths = copy;
        
        if (h1>1e7&&h2>1e7)continue;
        //int score=h1<h2?h1:h2;
        //int score=(h1+h2)/2;
        typedef long long LL;
        //double score=(min(h1,h2)+1)*1.0*200*200-minConf;
        double score=min(h1,h2);
        //if (useoracle==1 || useoracle==3)
        //score=h1<h2?h1:h2;
        if (useoracle==1)
        {
            if (maxScore < score)
            {
                chosenConflict=conflict;
                maxScore = score;
            }
        }else
        {
            if (minScore > score)
            {
                chosenConflict=conflict;
                minScore = score;
            }
        }
        if (maxScore < score)maxScore=score;
        if (minScore > score)minScore=score;
        
        
        
        
        if (collectingFeatures)
        {
            
            int K=5;
            //featuresCollected++;
            
            std::vector<double> features;
            features.clear();
            
            //score
            features.push_back(score);
            //features.push_back(min(h1,h2));

            //current makespan (1)
            //features.push_back(currentmakespan);//1
            //current depth
            //features.push_back(currentdepth);//2
            //current f,g,h value (3)
            //features.push_back(curr_fval);//3
            //features.push_back(curr_gval);//4
            //features.push_back(curr_hval);//5
            
            //number of conflicts(1)
            //features.push_back(node.num_of_collisions);//6

            
            
            
            //conflict type (2)
            if (get<3>(*node.conflict) < 0) //vertex conflict
                features.push_back(1),features.push_back(0);//1,2
            else features.push_back(0),features.push_back(1);
            
            //cardinality of the conflict (3)

            int conflictType=node.confType[numC];

            int conflicta=0,conflictb=0,conflictc=0;
            if(conflictType==1)conflicta=1;
            else if (conflictType==0)conflictb=1;
            else conflictc=1;
            features.push_back(conflicta);//3
            features.push_back(conflictb);//4
            features.push_back(conflictc);//5
            
            double pathL0=paths[get<0>(*node.conflict)]->size();
            double pathL1=paths[get<1>(*node.conflict)]->size();
            
            int u=get<2>(*node.conflict);
            int v=get<3>(*node.conflict);
            if (v==-1)v=u;
            
            //number of conflicts i,j resolved so far (2)
            features.push_back(countResolveforAgent[get<0>(*node.conflict)]); //6 x
            features.push_back(countResolveforAgent[get<1>(*node.conflict)]); //7 x
            
            features.push_back(countResolveforAgent[get<0>(*node.conflict)]+countResolveforAgent[get<1>(*node.conflict)]);//8
            
            //number of conflicts node u,v resolved so far (2)
            features.push_back(countResolveforNode[u]); //9
            features.push_back(countResolveforNode[v]); //10
            
            features.push_back(countResolveforNode[u]+countResolveforNode[v]);//11


            //number of conflicts i,j involved (2)
            features.push_back(conflictsInvolved[get<0>(*node.conflict)]); //12 x
            features.push_back(conflictsInvolved[get<1>(*node.conflict)]); //13 x
            features.push_back(conflictsInvolved[get<0>(*node.conflict)]+conflictsInvolved[get<1>(*node.conflict)]);//14

            //time of conflict (1)
            int timestep=get<4>(*node.conflict);
            features.push_back(timestep);//15
            
            //time of conflict/ currentmakespan (1)
            features.push_back(timestep*1.0/currentmakespan);//16
            
            //pathlength of agents i,j (2)
            features.push_back(pathL0);//17 ??
            features.push_back(pathL1);//18 ??
            features.push_back(pathL0+pathL1);//19
            features.push_back(abs(pathL0-pathL1));//20
            features.push_back(max(pathL0*1./pathL1,pathL1*1./pathL0));//21
            
            features.push_back(pathL0-shortestPathLen[get<0>(*node.conflict)]);//22
            features.push_back(pathL1-shortestPathLen[get<1>(*node.conflict)]);//23
            
            features.push_back(pathL0*1./shortestPathLen[get<0>(*node.conflict)]);//24
            features.push_back(pathL1*1./shortestPathLen[get<1>(*node.conflict)]);//25

            
            //pathlength / curr_fval (2)
            features.push_back(pathL0/curr_fval);//26
            features.push_back(pathL1/curr_fval);//27
            
            //whether i,j have reached goal (2)
            features.push_back(int(paths[get<0>(*node.conflict)]->size()>timestep));//28
            features.push_back(int(paths[get<1>(*node.conflict)]->size()>timestep));//29
            
            //pathlength-time (2)
            features.push_back(pathL0-timestep);//30
            features.push_back(pathL1-timestep);//31
            
            //(pathlength-time)/pathlength (2)
            features.push_back((pathL0-timestep)/pathL0);//32
            features.push_back((pathL1-timestep)/pathL1);//33
            
            //cerr<<"ok1"<<endl;
            
            //
            std::vector<int> nodeKstep(K,0);// nodes within k steps in the graph
            
            std::vector<int> distu;
            ml->findKNeighbor(u,v,K,distu,nodeKstep);
            //std::vector<int> distv=v==-1?distu:findKNeighbor(v,K);
            
            //cerr<<"ok2"<<endl;

            
            std::vector<int> conflictKstep(K+1,0);//conflicts within k steps in time-expanded graph
            std::vector<int> conflictKstep2(K+1,0);//conflicts within k steps in the graph
            std::vector<int> agentKstep(K+1,0);//agents within k steps in time-expanded graph
            
            for (std::shared_ptr<Conflict> tpCon:node.allConf)
            {
                int tu=get<2>(*tpCon);
                int tv=get<3>(*tpCon);
                int tt=get<4>(*tpCon);
                if (tv==-1)tv=tu;
                int tdist=min(distu[tu],distu[tv]);
                //int tdistv=min(distv[tu],distv[tv]);
                //int tdist=min(tdistu,tdistv);
                if(tdist>=1e7)continue;
                if(tdist+abs(tt-timestep)<=K)
                    conflictKstep[tdist+abs(tt-get<4>(*node.conflict))]+=1;
                if(tdist<=K)conflictKstep2[tdist]+=1;
            }
            
            //cerr<<"ok3"<<endl;

            for(int a=0;a<num_of_agents;a++)
            {
                //int loc_a=paths[a]->at(timestep).location;
                int loc_a;
                int len_a=paths[a]->size();
                //if (distu[loc_a]>1e7 && distv[loc_a]>1e7) continue;
                for(int tt=timestep;tt>=0&&tt>=timestep-K;tt--)
                {
                    if (tt<len_a) loc_a=paths[a]->at(tt).location;
                    else loc_a=paths[a]->at(len_a-1).location;
                        
                    if(timestep-tt+distu[loc_a]<=K)
                        agentKstep[timestep-tt+distu[loc_a]]+=1;
                    else break;
                }
                for(int tt=timestep+1;tt<=timestep+K;tt++)
                {
                    if (tt<len_a) loc_a=paths[a]->at(tt).location;
                    else loc_a=paths[a]->at(len_a-1).location;
                    if(tt-timestep+distu[loc_a]<=K)
                        agentKstep[tt-timestep+distu[loc_a]]+=1;
                    else break;
                }
            }
            //cerr<<"ok4"<<endl;
            
            
            // (4k+3)
            features.insert(features.end(), conflictKstep.begin(), conflictKstep.end()); //34-39
            features.insert(features.end(), conflictKstep2.begin(), conflictKstep2.end());//40-45
            features.insert(features.end(), agentKstep.begin(), agentKstep.end());//46-51
            features.insert(features.end(), nodeKstep.begin(), nodeKstep.end());//52-56
            
            int a1 = min(get<0>(*node.conflict), get<1>(*node.conflict));
            int a2 = max(get<0>(*node.conflict), get<1>(*node.conflict));
            for(int tt=timestep-2;tt<=timestep+2;tt++)
            {
                int mddw1=0,mddw2=0;
                if(tt>=0)
                {
                    if(tt<paths[a1]->size())mddw1=paths[a1]->at(tt).width;
                    if(tt<paths[a2]->size())mddw2=paths[a2]->at(tt).width;
                }
                features.push_back(mddw1);
                features.push_back(mddw2);
            }//3rd new add
            
            int idx = a1 * num_of_agents + a2;
            auto got = node.conflictGraph.find(idx);
            if (got != node.conflictGraph.end())
                features.push_back(got->second);
            else {features.push_back(0);cerr<<"no weight in DG!"<<endl;exit(0);}
            
            featureSet.push_back(features);
            featureScore.push_back(score);

            /*
            freopen(featureFile.c_str(),"a",stdout);
            for(auto feat:features)printf("%.4lf ",feat);
            puts("");
            fclose(stdout);*/
        }
        
        //releaseMDDMemory(curr->agent_id);
    }
    if (fabs(minScore-maxScore)<1e-7)sameScore+=1;
    else
    {
        freopen(featureFile.c_str(),"a",stdout);
        sort(featureScore.begin(),featureScore.end());
        if(useoracle==1)
            reverse(featureScore.begin(),featureScore.end());
        int tp=featureScore.size()/10;
        int tpj=tp+1;
        for(;tpj<featureScore.size();tpj++)
            if(featureScore[tp]!=featureScore[tpj])break;
        //printf("%d %d ",tp,tpj);
        if(tpj*1.0/featureScore.size()>=0.2)
        {
            int tpk=tp;
            for(;tp>0;tp--)
                if(featureScore[tpk]!=featureScore[tp])break;
            //printf("%d\n",tp);
        }
        //printf("\n");
        diffScore+=1;
        featuresCollected++;
        /*
        vector<double> norm(featureSize,0);
        for(auto features:featureSet)
        {
            for(int i=0;i<features.size();i++)
                norm[i]=max(norm[i],abs(features[i]));
        }
        for(int i=0;i<norm.size();i++)
            if(norm[i]<1e-8)norm[i]=1;
        for(auto features:featureSet)
        {
            printf("0 qid:1 ");
            for(int i=1;i<=features.size();i++)
                printf("%d:%.4lf ",i,features[i-1]/norm[i-1]);
            puts("");
        }*/
        for(auto feature:featureSet)
        {
            feature.insert(feature.begin(),featuresCollected);
            if (useoracle==1)
            {
                if(feature[1]+1e-8>=featureScore[tp])feature.insert(feature.begin(),1);
                else feature.insert(feature.begin(),0);
            } else
            {
                if(feature[1]-1e-8<=featureScore[tp])feature.insert(feature.begin(),1);
                else feature.insert(feature.begin(),0);
            }
            
            for(auto feat:feature)
                if(fabs(feat-int(feat+0.01))<1e-8)printf("%d ",int(feat+0.01));
                else printf("%.4lf ",feat);
            puts("");
            //for(int i=)
        }
        puts("");
        fclose(stdout);
    }
    node.conflict = chosenConflict;
    chooseConflictTime+=(std::clock() - t1) * 1000.0/ CLOCKS_PER_SEC;
}

void ICBSSearch::releaseMDDMemory(int id)
{
    if (id < 0 || mddTable.empty() || (int)mddTable[id].size() < max_num_of_mdds)
        return;
    int minLength = INT_MAX;
    for (auto mdd : mddTable[id])
    {
        if ((int)mdd.second->levels.size() < minLength)
            minLength = mdd.second->levels.size();
    }
    for (MDDTable::iterator mdd = mddTable[id].begin(); mdd != mddTable[id].end();)
    {
        if (mdd->second->levels.size() == minLength)
        {
            delete mdd->second;
            mdd = mddTable[id].erase(mdd);
            num_released_mdds++;
        }
        else
        {
            mdd++;
        }
    }
    return;
}


ICBSSearch::ICBSSearch(const MapLoader* ml, vector<SingleAgentICBS*>& search_engines, const vector<list<Constraint>>& constraints,
    vector<vector<PathEntry>>& paths_found_initially, double f_w, int initial_h,
    heuristics_type h_type, bool PC, bool rectangleReasoning,
    int cost_upperbound, double time_limit, int screen):
    focal_w(f_w), time_limit(time_limit), h_type(h_type), PC(PC), screen(screen), cost_upperbound(cost_upperbound),
    rectangleReasoning(rectangleReasoning), ml(ml),
    search_engines(search_engines), initial_constraints(constraints), paths_found_initially(paths_found_initially)
{
    
    HL_num_expanded = 0;
    HL_num_generated = 0;
    LL_num_expanded = 0;
    LL_num_generated = 0;
    
    num_of_agents = search_engines.size();
    
    countResolveforAgent.resize(num_of_agents,0);
    int rr=ml->rows;
    int cc=ml->cols;
    countResolveforNode.resize(rr*cc,0);
    
   
    
    solution_found = false;
    solution_cost = -2;

    // generate dummy start and update data structures
    dummy_start = new ICBSNode();
    dummy_start->agent_id = -1;
    dummy_start->g_val = 0;
    paths.resize(num_of_agents);
    for (int i = 0; i < num_of_agents; i++)
    {
        paths[i] = &paths_found_initially[i];
        dummy_start->makespan = max(dummy_start->makespan, paths_found_initially[i].size() - 1);
        dummy_start->g_val += paths_found_initially[i].size() - 1;
    }
    dummy_start->h_val = initial_h;
    dummy_start->f_val = dummy_start->g_val;

    dummy_start->depth = 0;

    dummy_start->open_handle = open_list.push(dummy_start);
    dummy_start->focal_handle = focal_list.push(dummy_start);

    HL_num_generated++;
    dummy_start->time_generated = HL_num_generated;
    allNodes_table.push_back(dummy_start);
    findConflicts(*dummy_start);

    min_f_val = dummy_start->f_val;
    focal_list_threshold = min_f_val * focal_w;

    if (rectangleReasoning)
        mddTable.resize(num_of_agents);
}

ICBSSearch::ICBSSearch(const MapLoader& ml, const AgentsLoader& al, double f_w, heuristics_type h_type,
    bool PC, bool rectangleReasoning,
    double time_limit, int screen):
    focal_w(f_w), time_limit(time_limit), h_type(h_type), PC(PC), screen(screen),
    rectangleReasoning(rectangleReasoning), ml(&ml),
    num_of_agents(al.num_of_agents)
{

    countResolveforAgent.resize(num_of_agents,0);
    int rr=ml.rows;
    int cc=ml.cols;
    countResolveforNode.resize(rr*cc,0);
    
    initial_constraints.resize(num_of_agents);

    search_engines = vector < SingleAgentICBS* >(num_of_agents);
    for (int i = 0; i < num_of_agents; i++)
    {
        int init_loc = ml.linearize_coordinate((al.initial_locations[i]).first, (al.initial_locations[i]).second);
        int goal_loc = ml.linearize_coordinate((al.goal_locations[i]).first, (al.goal_locations[i]).second);
        ComputeHeuristic ch(init_loc, goal_loc, ml.get_map(), ml.rows, ml.cols, ml.get_moves_offset());
        search_engines[i] = new SingleAgentICBS(init_loc, goal_loc, ml.get_map(), ml.rows*ml.cols,
            ml.get_moves_offset(), ml.cols);
        ch.getHVals(search_engines[i]->my_heuristic);
    }

    dummy_start = new ICBSNode();
    dummy_start->agent_id = -1;
    //cerr<<"okk"<<endl;

    
    // initialize paths_found_initially
    paths.resize(num_of_agents, NULL);
    shortestPathLen.resize(num_of_agents,0);
    paths_found_initially.resize(num_of_agents);
    vector < list< pair<int, int> > > cons_vec;
    for (int i = 0; i < num_of_agents; i++)
    {
        CAT cat(dummy_start->makespan + 1);  // initialized to false
        updateReservationTable(cat, i, *dummy_start);

        if (search_engines[i]->findPath(paths_found_initially[i], cons_vec, cat, 0) == false)
            cout << "NO SOLUTION EXISTS";

        paths[i] = &paths_found_initially[i];
        shortestPathLen[i]=paths[i]->size();

        dummy_start->makespan = max(dummy_start->makespan, paths_found_initially[i].size() - 1);

        LL_num_expanded += search_engines[i]->num_expanded;
        LL_num_generated += search_engines[i]->num_generated;
    }
    //cerr<<"okk"<<endl;



    // generate dummy start and update data structures
    dummy_start->g_val = 0;
    for (int i = 0; i < num_of_agents; i++)
        dummy_start->g_val += paths[i]->size() - 1;
    dummy_start->h_val = 0;
    dummy_start->f_val = dummy_start->g_val;

    dummy_start->depth = 0;
    
    dummy_start->open_handle = open_list.push(dummy_start);
    dummy_start->focal_handle = focal_list.push(dummy_start);

    HL_num_generated++;
    dummy_start->time_generated = HL_num_generated;
    allNodes_table.push_back(dummy_start);
    findConflicts(*dummy_start);

    min_f_val = dummy_start->f_val;
    focal_list_threshold = min_f_val * focal_w;
    //if(h_type == heuristics_type::DG || h_type == heuristics_type::PAIR)
    //    dummy_start->conflictGraph.resize(num_of_agents * num_of_agents, -1);
    /*if (h_type == heuristics_type::DG && !EPEA4PAIR)
        mdds_initially.resize(num_of_agents);*/

    if (screen >= 2) // print start and goals
    {
        al.printAgentsInitGoal();
    }

    hTable.resize(num_of_agents);
    for (int i = 0; i < num_of_agents; i++)
    {
        hTable[i].resize(num_of_agents);
    }

    if (rectangleReasoning || h_type == heuristics_type::DG || h_type == heuristics_type::WDG)
        mddTable.resize(num_of_agents);
}

inline void ICBSSearch::releaseClosedListNodes()
{
    for (list<ICBSNode*>::iterator it = allNodes_table.begin(); it != allNodes_table.end(); it++)
        delete *it;
}

inline void ICBSSearch::releaseOpenListNodes()
{
    while(!open_list.empty())
    {
        ICBSNode* curr = open_list.top();
        open_list.pop();
        delete curr;
    }
}

ICBSSearch::~ICBSSearch()
{
    releaseClosedListNodes();
    /*if (!mdds_initially.empty())
    {
        for (auto mdd : mdds_initially)
        {
            if (mdd != NULL)
                delete mdd;
        }
    }*/
    if (!mddTable.empty())
    {
        for (int i = 0; i < num_of_agents; i++)
        {
            for (auto mdd : mddTable[i])
            {
                    delete mdd.second;
            }
        }
    }
}

void ICBSSearch::clearSearchEngines()
{
    for (size_t i = 0; i < search_engines.size(); i++)
        delete (search_engines[i]);
}


bool ICBSSearch::validateSolution() const
{
    for (int a1 = 0; a1 < num_of_agents; a1++)
    {
        for (int a2 = a1 + 1; a2 < num_of_agents; a2++)
        {
            size_t min_path_length = paths[a1]->size() < paths[a2]->size() ? paths[a1]->size() : paths[a2]->size();
            for (size_t timestep = 0; timestep < min_path_length; timestep++)
            {
                int loc1 = paths[a1]->at(timestep).location;
                int loc2 = paths[a2]->at(timestep).location;
                if (loc1 == loc2)
                {
                    std::cout << "Agents "  << a1 << " and " << a2 << " collides at " << loc1 << " at timestep " << timestep << std::endl;
                    return false;
                }
                else if (timestep < min_path_length - 1
                    && loc1 == paths[a2]->at(timestep + 1).location
                    && loc2 == paths[a1]->at(timestep + 1).location)
                {
                    std::cout << "Agents " << a1 << " and " << a2 << " collides at (" <<
                        loc1 << "-->" << loc2 << ") at timestep " << timestep << std::endl;
                    return false;
                }
            }
            if (paths[a1]->size() != paths[a2]->size())
            {
                int a1_ = paths[a1]->size() < paths[a2]->size() ? a1 : a2;
                int a2_ = paths[a1]->size() < paths[a2]->size() ? a2 : a1;
                int loc1 = paths[a1_]->back().location;
                for (size_t timestep = min_path_length; timestep < paths[a2_]->size(); timestep++)
                {
                    int loc2 = paths[a2_]->at(timestep).location;
                    if (loc1 == loc2)
                    {
                        std::cout << "Agents " << a1 << " and " << a2 << " collides at " << loc1 << " at timestep " << timestep << std::endl;
                        return false; // It's at least a semi conflict
                    }
                }
            }
        }
    }
    return true;
}
