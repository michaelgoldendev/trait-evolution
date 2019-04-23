using FastaIO
using Nullables

export NodeData
abstract type NodeData end

export EmptyNodeData
mutable struct EmptyNodeData <: NodeData

end

export MyNodeData
mutable struct MyNodeData <: NodeData
  sequence::AbstractString

  MyNodeData(sequence::AbstractString) = new(sequence)
end

export TreeNode
mutable struct TreeNode
  parent::Nullable{TreeNode}
  children::Array{TreeNode,1}
  branchlength::Float64
  name::AbstractString
  data::NodeData
  nodeindex::Int
  seqindex::Int

  TreeNode(branchlength::Float64, name::AbstractString) = new(Nullable{TreeNode}(),TreeNode[],branchlength,name, EmptyNodeData(),0,0)
  TreeNode(branchlength::Float64, name::AbstractString, data::NodeData) = new(Nullable{TreeNode}(),TreeNode[],branchlength,name,data,0,0)
end

Base.first(node::TreeNode) = 1
Base.isdone(node::TreeNode,state) = length(node.children) == state-1
function Base.iterate(node::TreeNode,state::Int=1)
    if state <= length(node.children)
        return node.children[state], state+1
    else
        return nothing
    end
end
Base.length(node::TreeNode) = length(node.children)

export roottree
function roottree(root::TreeNode,index::Int=1)
  newroot = TreeNode(0.0, "")
  child = splice!(root.children,index)
  child.parent = Nullable{TreeNode}()
  addchild(newroot,child)
  br = child.branchlength/2.0
  child.branchlength = br
  root.branchlength = br
  addchild(newroot,root)
  return newroot
end

export isroot
function isroot(node::TreeNode)
  return isnull(node.parent)
end

export isleafnode
function isleafnode(node::TreeNode)
  return length(node.children) == 0
end

export addchild
function addchild(parent::TreeNode, child::TreeNode)
  if isnull(child.parent)
    push!(parent.children, child)
    child.parent = Nullable{TreeNode}(parent)
  else
    if insubtree(parent, child)
      throw(ArgumentError("Cannot move node to a subtree of itself!"))
    else
      children = get(child.parent).children
      index = findfirst(children, child)
      deleteat!(children, index)
      child.parent = Nullable{TreeNode}(parent)
      push!(parent.children, child)
    end
  end
end

export insubtree
function insubtree(node::TreeNode, subtree::TreeNode)
  if subtree == node
    return true
  else
    for child in subtree
      if insubtree(node, child)
        return true
      end
    end
  end
  return false
end


function getnewickhelper(node::TreeNode)
  if length(node.children) == 0
    return string(node.name,":", node.branchlength)
  else
    ret = join(AbstractString[getnewickhelper(child) for child in node], ",")
    return string("(",ret,")",node.name,":", node.branchlength)
  end
end

export getnewick
function getnewick(node::TreeNode)
  return string(getnewickhelper(node),";")
end

function treefromnewickhelper(newick::AbstractString)
  startindex = 1
  endindex = length(newick)
  lm = match(r"[\)][^\)]*$", newick)

  tag = newick
  childstrings = AbstractString[]
  if lm != nothing
    tag = newick[lm.offset+1:end]
    childstring = newick[startindex+1:lm.offset-1]

    child = ""
    a = 0
    for i=1:length(childstring)
      if childstring[i] == '('
        a += 1
        child = string(child,childstring[i])
      elseif childstring[i] == ')'
        a -= 1
        child = string(child,childstring[i])
      elseif childstring[i] == ',' && a == 0
        push!(childstrings, child)
        child = ""
      else
        child = string(child,childstring[i])
      end
    end
    if child != ""
      push!(childstrings, child)
    end
  end
  spl = split(tag,":")
  name = ""
  branchlength = 0.0
  if length(spl) > 0
    name = strip(spl[1])
  end
  if length(spl) > 1
    branchlength = parse(Float64,spl[2])
  end

  return childstrings,name,branchlength
end

export gettreefromnewick
"""
  gettreefromnewick(newick)

Returns a TreeNode object from `newick` string.
"""
function gettreefromnewick(newick::AbstractString)
  childstrings,name,branchlength = treefromnewickhelper(rstrip(strip(newick),';'))
  node = TreeNode(branchlength,name, EmptyNodeData())
  for childstring in childstrings
    addchild(node, gettreefromnewick(childstring))
  end
  return node
end

export prettyprintstring
function prettyprintstring(node::TreeNode, spaces::Int=0)
  ret = string(repeat("----",spaces),"+++++ ", node.name, "(", node.branchlength,")", "\n")
  for child in node
      ret = string(ret, prettyprintstring(child,spaces+1))
  end
  return ret
end

export ExtraNodeData
mutable struct ExtraNodeData <: NodeData
  nodeindex::Int
  seqindex::Int
  ExtraNodeData(seqindex::Int) = new(0,seqindex)
end

export annotatetree
function annotatetree(node::TreeNode, seqnametoindex::Dict{AbstractString,Int}, nodelist::Array{TreeNode,1}=TreeNode[])
  push!(nodelist,node)
  node.data = ExtraNodeData(get(seqnametoindex,node.name,0))
  node.data.nodeindex = length(nodelist)
  node.nodeindex = length(nodelist)
  node.seqindex = get(seqnametoindex,node.name,0)
  for childnode in node
    annotatetree(childnode,seqnametoindex,nodelist)
  end
end

export getnodelist
function getnodelist(node::TreeNode, nodelist::Array{TreeNode,1}=TreeNode[])
  push!(nodelist,node)
  for childnode in node
    getnodelist(childnode,nodelist)
  end
  return nodelist
end

export getpattern
function getpattern(data::Array{Float64,3}, node::TreeNode, col::Int, pattern::Array{Int8,1}=Int8[])
  if isleafnode(node)
    for b in data[node.seqindex,col,:]
      push!(pattern, Int8(b))
    end
  else
    for childnode in node
      getpattern(data,childnode,col,pattern)
    end
  end

  return (node.nodeindex, pattern)
end

function countnodes(node::TreeNode)
  if isleafnode(node)
    return 1
  else
    s = 0
    for chilnode in node
      s += countnode(childnode)
    end
    return 1+s
  end
end

export getnodedepth
function getnodedepth(node::TreeNode)
  depth = 1
  current = node
  while !isnull(current.parent)
    current = get(current.parent)
    depth += 1
  end
  return depth
end

export rootnodedistance
function rootnodedistance(node::TreeNode)
  depth = 0.0
  current = node
  while !isnull(current.parent)
    depth += current.branchlength
    current = get(current.parent)
  end
  return depth
end

export countleafnodes
function countleafnodes(node::TreeNode)
  if isleafnode(node)
    return 1
  else
    s = 0
    for childnode in node
      s += countleafnodes(childnode)
    end
    return s
  end
end
