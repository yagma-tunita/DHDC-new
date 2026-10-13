#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <queue>
#include <map>

struct node
{
	int dimension;
	int pathA;	
	int pathB;	
	int tag;
	bool tail;
	// std::vector<bool> pathA;	// A
	// std::vector<bool> pathB;	// B

	node(int dim, int A, int B, int t, bool i)
		: dimension(dim), pathA(A), pathB(B), tag(t), tail(i) {
	}
	// node(int dim, const std::vector<bool>& A, const std::vector<bool>& B, int t, bool i)
	// 	: dimension(dim), pathA(A), pathB(B), tag(t), tail(i) {
	// }

	bool operator==(const node& other) const
	{
		return dimension == other.dimension
			&& pathA == other.pathA
			&& pathB == other.pathB
			&& tail == other.tail
			&& tag == other.tag;
	}
	bool operator!=(const node& other) const { return !(*this == other); }
	bool operator<(const node& other)const
	{
		if (this->dimension != other.dimension)
			return this->dimension < other.dimension;
		if (this->pathA != other.pathA)
			return this->pathA < other.pathA;
		if (this->pathB != other.pathB)
			return this->pathB < other.pathB;
		if (this->tail != other.tail)
			return this->tail < other.tail;
		if (this->tag != other.tag)
			return this->tag < other.tag;
		return false;
	}
	struct ncmp
	{
		bool operator()(const node& a, const node& b) const
		{
			return a < b;
		}
	};

	using dhmap = std::map<node, node, ncmp>;

	void print_node()
	{
		std::cout << '(';
		for (int i = 0; i < this->dimension - 1; i++)
			std::cout << ((this->pathA >> i) & 1);
		for (int i = 0; i < this->dimension - 1; i++)
			std::cout << ((this->pathB >> i) & 1);
		if (this->tail)
			std::cout << 1;
		else
			std::cout << 0;
		std::cout << ',' << this->tag << ')';
	}
};

struct NodeHash
{
	std::size_t operator()(const node& n) const
	{
		std::size_t h = std::hash<int>()(n.dimension);
		h ^= std::hash<int>()(n.pathA) + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<int>()(n.pathB) + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<int>()(n.tag) + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<bool>()(n.tail) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};

struct HsdicNode
{
	int dimension;
	std::vector<bool> bits;
	int tag;
	HsdicNode() {}
	HsdicNode(int d, const std::vector<bool>& b, int y)
		: dimension(d), bits(b), tag(y) {}
	bool operator==(const HsdicNode& o) const
	{
		return dimension == o.dimension && bits == o.bits && tag == o.tag;
	}
	bool operator!=(const HsdicNode& o) const { return !(*this == o); }
};

struct HsdicNodeHash
{
	std::size_t operator()(const HsdicNode& n) const
	{
		std::size_t h = std::hash<int>()(n.dimension);
		h ^= std::hash<int>()(n.tag) + 0x9e3779b9 + (h << 6) + (h >> 2);
		std::size_t hb = 0;
		for (bool b : n.bits) hb = hb * 131 + (b ? 1 : 0);
		h ^= hb + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};

// BCube
struct BcubeNode
{
	int n;
	int k;
	std::vector<int> addr;
	BcubeNode() {}
	BcubeNode(int nn, int kk, const std::vector<int>& a)
		: n(nn), k(kk), addr(a) {}
	bool operator==(const BcubeNode& o) const
	{
		return n == o.n && k == o.k && addr == o.addr;
	}
	bool operator!=(const BcubeNode& o) const { return !(*this == o); }
};

struct BcubeNodeHash
{
	std::size_t operator()(const BcubeNode& x) const
	{
		std::size_t h = std::hash<int>()(x.n);
		h ^= std::hash<int>()(x.k) + 0x9e3779b9 + (h << 6) + (h >> 2);
		for (int v : x.addr) h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};

// DCell(n,k)
struct DcellNode
{
	int n;
	int k;
	std::vector<int> addr;
	DcellNode() {}
	DcellNode(int nn, int kk, const std::vector<int>& a)
		: n(nn), k(kk), addr(a) {}
	bool operator==(const DcellNode& o) const
	{
		return n == o.n && k == o.k && addr == o.addr;
	}
	bool operator!=(const DcellNode& o) const { return !(*this == o); }
};

struct DcellNodeHash
{
	std::size_t operator()(const DcellNode& x) const
	{
		std::size_t h = std::hash<int>()(x.n);
		h ^= std::hash<int>()(x.k) + 0x9e3779b9 + (h << 6) + (h >> 2);
		for (int v : x.addr) h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};

// Fat-Tree(k)
struct FatTreeNode
{
	int k;
	int pod;
	int edgeSwitch;
	int server;
	FatTreeNode() {}
	FatTreeNode(int kk, int p, int e, int s)
		: k(kk), pod(p), edgeSwitch(e), server(s) {}
	bool operator==(const FatTreeNode& o) const
	{
		return k == o.k && pod == o.pod && edgeSwitch == o.edgeSwitch && server == o.server;
	}
	bool operator!=(const FatTreeNode& o) const { return !(*this == o); }
};

struct FatTreeNodeHash
{
	std::size_t operator()(const FatTreeNode& x) const
	{
		std::size_t h = std::hash<int>()(x.k);
		h ^= std::hash<int>()(x.pod) + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<int>()(x.edgeSwitch) + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<int>()(x.server) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};

// FiConn(k)
struct FiconNode
{
	int k;
	std::vector<int> addr;
	FiconNode() {}
	FiconNode(int kk, const std::vector<int>& a)
		: k(kk), addr(a) {}
	bool operator==(const FiconNode& o) const
	{
		return k == o.k && addr == o.addr;
	}
	bool operator!=(const FiconNode& o) const { return !(*this == o); }
};

struct FiconNodeHash
{
	std::size_t operator()(const FiconNode& x) const
	{
		std::size_t h = std::hash<int>()(x.k);
		for (int v : x.addr) h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};
