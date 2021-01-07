#include <iostream>
#include <tuple>

// Given a binary tree, process a number of queries with parameters: A, d.
// For each query, answer with the value of node B if there is such a node B
// that its distance from A equals d, and answer with -1 otherwise.
// Preprocessing the tree is allowed (it is given in the form of a sequence
// of parent-child relationships).

// The solution consists in memorizing, for each node A, the node B that is farthest
// from it, and the distance from A to B. Then for each query, if d is smaller than
// that distance, the appropriate node is looked for on the route from A to B
// (which is simplified by first finding their lowest common ancestor). 
// this results in logarithmic time of answering queries.

// ceil of log n
int log_ceil(int n){
    int res = 0;
    int pow = 1;
    while(pow < n){
        pow *= 2;
        res++;
    }
    return res;
}

class node{
    int level;
    int id;
    node* left_son;
    node* right_son;
    node** ancestors;

    node* farthest_left_leaf;
    int farthest_left_leaf_dist;

    node* farthest_right_leaf;
    int farthest_right_leaf_dist;

    node* farthest;
    int farthest_dist;

public:
    int get_id(){return id;}
    node* get_lson(){return left_son;}
    node* get_rson(){return right_son;}
    int get_level(){return level;}
    int get_farthest_dist(){return farthest_dist;}
    node* get_farthest(){return farthest;}
    node** get_ancestors(){return ancestors;}
    static int max_height;
    node(int id){
        this->id = id;
        this->ancestors = (node**) malloc(max_height * sizeof(node*));
        for(int i = 0; i < max_height; i++) ancestors[i] = nullptr;
    }
    void set_lson(int n, node* nodes){this->left_son = n == -1 ? nullptr : &nodes[n-1];}
    void set_rson(int n, node* nodes){this->right_son = n == -1 ? nullptr : &nodes[n-1];}
    void set_father(int n, node* nodes){
        this->ancestors[0] = &nodes[n];
    }
    void set_level(int l){this->level = l;}
    void fill_aux(int level);
    void fill_leaf_dists();
    void set_farthest_left_leaf(node* lleaf){this->farthest_left_leaf = lleaf;}
    void set_farthest_left_leaf_dist(int ldist){this->farthest_left_leaf_dist = ldist;}
    void set_farthest_right_leaf(node* rleaf){this->farthest_right_leaf = rleaf;}
    void set_farthest_right_leaf_dist(int rdist){this->farthest_right_leaf_dist = rdist;}
    void print_node(int m);
    void fill_cousin_dists_aux(int f, node* fathers_farthest);
    void fill_cousin_dists();
};

int node::max_height;

void node::fill_aux(int level){
    int i = 1;
    this->level = level;
    while(ancestors[i-1]->ancestors[i-1] != nullptr){
        ancestors[i] = ancestors[i-1]->ancestors[i-1];
        i++;
    }
    if(this->get_lson() != nullptr){
        this->get_lson()->fill_aux(level+1);
    }
    if(this->get_rson() != nullptr){
        this->get_rson()->fill_aux(level+1);
    }
}

void fill_ancestry_and_stuff(node* tr){
    tr->set_level(0);
    if(tr->get_lson() != nullptr) tr->get_lson()->fill_aux(1);
    if(tr->get_rson() != nullptr) tr->get_rson()->fill_aux(1);
}

node* furthest_leaf(node* self, int* dist){
    if(self == nullptr) return nullptr;
    else if(self->get_lson() == nullptr && self->get_rson() == nullptr){
        self->set_farthest_left_leaf_dist(0);
        self->set_farthest_right_leaf_dist(0);
        self->set_farthest_left_leaf(nullptr);
        self->set_farthest_right_leaf(nullptr);
        *dist = 0;
        return self;
    }
    int dl = 0;
    node* left_leaf = furthest_leaf(self->get_lson(), &dl);
    self->set_farthest_left_leaf(left_leaf);
    if(left_leaf != nullptr){
        dl++;
        self->set_farthest_left_leaf_dist(dl);
    }
    else self->set_farthest_left_leaf_dist(-1);
    int dr = 0;
    node* right_leaf = furthest_leaf(self->get_rson(), &dr);
    self->set_farthest_right_leaf(right_leaf);
    if(right_leaf != nullptr){
        dr++;
        self->set_farthest_right_leaf_dist(dr);
    }
    else self->set_farthest_right_leaf_dist(-1);
    if(dl > dr){
        *dist = dl;
        return left_leaf;
    }
    else{
        *dist = dr;
        return right_leaf;
    }

}

void node::fill_leaf_dists(){
    int dist = 0;
    this->farthest_left_leaf = furthest_leaf(this->get_lson(), &dist);
    this->farthest_left_leaf_dist = dist + 1;
    dist = 0;
    this->farthest_right_leaf = furthest_leaf(this->get_rson(), &dist);
    this->farthest_right_leaf_dist = dist + 1;
}


void node::fill_cousin_dists_aux(int f, node* fathers_farthest) {
    if(f >= farthest_right_leaf_dist && f >= farthest_left_leaf_dist){
        farthest = fathers_farthest;
        farthest_dist = f;
    }
    else if(farthest_right_leaf_dist >= f && farthest_right_leaf_dist >= farthest_left_leaf_dist){
        farthest = farthest_right_leaf;
        farthest_dist = farthest_right_leaf_dist;
    }
    else if(farthest_left_leaf_dist >= f && farthest_left_leaf_dist >= farthest_right_leaf_dist){
        farthest = farthest_left_leaf;
        farthest_dist = farthest_left_leaf_dist;
    }
    if(left_son != nullptr){
        if(f >= farthest_right_leaf_dist) left_son->fill_cousin_dists_aux(f+1, fathers_farthest);
        else left_son->fill_cousin_dists_aux(farthest_right_leaf_dist+1, farthest_right_leaf);
    }
    if(right_son != nullptr){
        if(f >= farthest_left_leaf_dist) right_son->fill_cousin_dists_aux(f+1, fathers_farthest);
        else right_son->fill_cousin_dists_aux(farthest_left_leaf_dist+1, farthest_left_leaf);
    }
}

void node::fill_cousin_dists(){
    node* left = get_lson();
    node* right = get_rson();
    if(farthest_left_leaf_dist > farthest_right_leaf_dist){
        farthest = farthest_left_leaf;
        farthest_dist = farthest_left_leaf_dist;
    }
    else if(farthest_right_leaf_dist != -1){
        farthest = farthest_right_leaf;
        farthest_dist = farthest_right_leaf_dist;
    }
    else{
        farthest = nullptr;
        farthest_dist = 0;
    }
    if(right == nullptr){
        if(left != nullptr) left->fill_cousin_dists_aux(1, this);
    }
    else {
        if (left != nullptr)
            left->fill_cousin_dists_aux(farthest_right_leaf_dist + 1, farthest_right_leaf);
    }
    if(left == nullptr){
        if(right != nullptr) right->fill_cousin_dists_aux(1, this);
    }
    else {
        if (right != nullptr)
            right->fill_cousin_dists_aux(farthest_left_leaf_dist + 1, farthest_left_leaf);
    }
}

node* jump_upwards(node* a, int n){
    if(n > a->get_level()) return nullptr;
    int i = 0;
    node* c = a;
    while(n != 0){
        if(n % 2 == 1) c = c->get_ancestors()[i];
        n /= 2;
        i++;
    }
    return c;
}

int pow(int a, int b){
    int res = 1;
    for(int i = 0; i < b; i++) res *= a;
    return res;
}

void lca(node* a, node* b, int* jumped_a, int* jumped_b){
    int a_level = a->get_level();
    int b_level = b->get_level();
    if(a_level < b_level){
        b = jump_upwards(b, b_level - a_level);
        *jumped_b += b_level - a_level;
    }
    else if(a_level > b_level){
        a = jump_upwards(a, a_level - b_level);
        *jumped_a += a_level - b_level;
    }
    if(b == a) return;
    int i = 0;
    while(a != b){
        i = 0;
        while(a->get_ancestors()[i] != nullptr && a->get_ancestors()[i] != b->get_ancestors()[i]) i++;
        if(i != 0)i -= 1;
        a = a->get_ancestors()[i];
        b = b->get_ancestors()[i];
        *jumped_a += pow(2, i);
        *jumped_b += pow(2, i);
    }
}

node* construct_tree(){
    int n;
    std::cin >> n;
    node::max_height = log_ceil(n);
    node* nodes = (node*) malloc(n * sizeof(node));
    for(int i = 0; i < n; i++){
        nodes[i] = node(i+1);
    }
    int left;
    int right;
    for(int i = 0; i < n; i++){
        std::cin >> left;
        std::cin >> right;
        nodes[i].set_lson(left, nodes);
        nodes[i].set_rson(right, nodes);
        if(left != -1) nodes[left-1].set_father(i, nodes);
        if(right != -1) nodes[right-1].set_father(i, nodes);
    }
    fill_ancestry_and_stuff(&nodes[0]);
    nodes[0].fill_leaf_dists();
    nodes[0].fill_cousin_dists();
    return nodes;
}

void process_query(node* tr){
    int no;
    std::cin >> no;
    int d;
    node* a = &tr[no-1];
    std::cin >> d;
    if(d == 0){
        std::cout << a->get_id() << "\n";
        return;
    }
    if(d > a->get_farthest_dist()){
        std::cout << -1 << "\n";
        return;
    }
    node* b = a->get_farthest();
    if(b == nullptr){
        std::cout << -1 << "\n";
        return;
    }
    int jumped_a = 0;
    int jumped_b = 0;
    lca(a, b, &jumped_a, &jumped_b);
    if(jumped_a >= d){
        a = jump_upwards(a, d);
        std::cout << a->get_id() << "\n";
    }
    else{
        b = jump_upwards(b, jumped_b - (d - jumped_a));
        std::cout << b->get_id() << "\n";
    }
}

void process_queries(node* tr){
    int m;
    std::cin >> m;
    for(int i = 0; i < m; i++) process_query(tr);
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    node* tr = construct_tree();
    process_queries(tr);
    return 0;
}