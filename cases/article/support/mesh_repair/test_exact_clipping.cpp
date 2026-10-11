#include "ArticleExactClippedVoronoiBuilder2D.hpp"
#include <VoronoiMeshMaker/Boundary2D/Builders/Boundary2DBuilder.hpp>
#include <VoronoiMeshMaker/Boundary2D/Shapes/Rectangle.hpp>
#include <VoronoiMeshMaker/Sites2D/Factory/SiteFactory.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiBuilder2D.hpp>
#include <iostream>
#include <iomanip>
int main(int argc,char** argv){try{auto n=argc>1?std::stoull(argv[1]):87552;auto b=vmm::b2d::make_boundary(vmm::b2d::Rectangle(vmm::b2d::Point2{-15,-5},30,10));auto s=vmm::s2d::make_sites(b,vmm::s2d::UniformRandom2D{n,12345},vmm::s2d::SiteValidationOptions{});auto d=vmm::vd2d::ArticleExactClippedVoronoiBuilder2D::build(s,b);long double area=0;std::size_t internal=0;for(auto& c:d.cells){area+=c.signed_area();for(auto& e:c.edges)if(!e.is_boundary_edge){++internal;if(e.neighbour_site_id.value<0)throw std::runtime_error("missing neighbour");}}if(n==100){
 auto original=vmm::vd2d::ClippedVoronoiBuilder2D::build(s,b);
 double max_area_difference=0;
 for(std::size_t i=0;i<d.cells.size();++i){
  const auto& x=d.cells[i];const auto& y=original.cells[i];
  max_area_difference=std::max(max_area_difference,std::abs(x.signed_area()-y.signed_area()));
  std::vector<long> a,z;
  for(const auto& e:x.edges)if(!e.is_boundary_edge)a.push_back(e.neighbour_site_id.value);
  for(const auto& e:y.edges)if(!e.is_boundary_edge)z.push_back(e.neighbour_site_id.value);
  std::sort(a.begin(),a.end());std::sort(z.begin(),z.end());
  if(a!=z)throw std::runtime_error("small mesh adjacency changed");
 }
 if(max_area_difference>1e-10)throw std::runtime_error("small mesh area changed");
 // Move one internal vertex in one cell only, while preserving endpoint/polygon
 // consistency. A real connectivity error must still be rejected.
 bool changed=false;
 for(auto& cell:d.cells){
  for(std::size_t e=0;e<cell.edges.size();++e)if(!cell.edges[e].is_boundary_edge){
   cell.polygon[e].x+=1e-5;cell.edges[e].a=cell.polygon[e];
   cell.edges[(e+cell.edges.size()-1)%cell.edges.size()].b=cell.polygon[e];changed=true;break;
  }
  if(changed)break;
 }
 bool rejected=false;
 try{vmm::vd2d::VoronoiFaceConnectivity2D::rebuild(d.cells,s,0);}catch(const std::runtime_error&){rejected=true;}
 if(!changed||!rejected)throw std::runtime_error("corrupted mesh was accepted");
 std::cout<<"REGRESSION adjacency unchanged; max_area_difference="<<max_area_difference<<"; corrupted mesh rejected\n";
 }
 std::cout<<std::setprecision(17)<<"ACCEPTED cells="<<d.cells.size()<<" internal_directed_faces="<<internal<<" area="<<area<<" error="<<std::abs(area-300)<<std::endl;if(std::abs(area-300)>1e-8)throw std::runtime_error("area mismatch");}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
