/* Copyright (C) 2010 Oliver Lemke <olemke@core-dump.info>
 
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any
 later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 USA. */

/*!
 \file   docserver.cc
 \author Oliver Lemke <olemke@core-dump.info>
 \date   2010-09-21
 
 \brief  Implementation of the arts documentation server.
 */

#include "docserver.h"

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
//#include <unistd.h>
#include "libmicrohttpd/platform.h"
#include "libmicrohttpd/microhttpd.h"
#include "messages.h"
#include "methods.h"
#include "workspace_ng.h"
#include "agenda_record.h"

#define DOCSERVER_NAME "ARTS built-in documentation server"

void limit_line_length (ostream& os,
                        ostringstream& curline,
                        ostringstream& token,
                        const String& indent,
                        size_t linelen);

void get_short_wsv_description(String &s, const String &desc);

bool format_paragraph (String &s, const String &indent, const size_t linelen,
                       const size_t offset = 0);

vector<string> &split(const string &s, char delim, vector<string> &elems) {
  stringstream ss(s);
  string item;
  while(getline(ss, item, delim))
  {
    elems.push_back(item);
  }
  return elems;
}

void ds_begin_page (ostream &os)
{
  os
  << "<!DOCTYPE html>" << endl
  << "<html>" << endl
  << "<head>" << endl
  << "<meta charset=\"utf-8\">" << endl
  << "<meta name=\"viewport\" content=\"width=device-width\" />" << endl
  << "<link rel=\"stylesheet\" href=\"/styles.css\">" << endl
  << "</head>" << endl
  << "<body>" << endl;
}

void ds_end_page (ostream &os)
{
  os << "</body>" << endl << "</html>";
}

String ds_insert_agenda_link (const String &aname)
{
  ostringstream ret;
  ret << "<a href=\"/agendas/" << aname << "\">" << aname << "</a>";
  return ret.str();
}

String ds_insert_group_link (const String &gname)
{
  ostringstream ret;
  ret << "<a href=\"/groups/" << gname << "\">" << gname << "</a>";
  return ret.str();
}

String ds_insert_wsm_link (const String &mname)
{
  ostringstream ret;
  ret << "<a href=\"/methods/" << mname << "\">" << mname << "</a>";
  return ret.str();
}

String ds_insert_wsv_link (const String &vname)
{
  ostringstream ret;
  
  // Find wsv id:
  map<String, Index>::const_iterator it = Workspace::WsvMap.find(vname);
  if ( it != Workspace::WsvMap.end() )
  {
    if (Workspace::wsv_data[it->second].Group() == get_wsv_group_id("Agenda"))
      ret << "<a href=\"/agendas/" << vname << "\">" << vname << "</a>";
    else
      ret << "<a href=\"/variables/" << vname << "\">" << vname << "</a>";
  }
  
  return ret.str();
}

void ds_list_agendas (ostream &os)
{
  Index i;
  
  os << "<tr><td colspan=\"2\"><a name=\"agendas\" />"
  << "<h2>Agendas</h2></td></tr>" << endl
  << "<tr><td><ul>" << endl;
  
  Index hitcount = 0;
  for ( i=0; i<Workspace::wsv_data.nelem(); ++i )
  {
    if (Workspace::wsv_data[i].Group() == get_wsv_group_id("Agenda"))
      hitcount++;
  }
  
  Index hitcount2 = 0;
  for ( i=0; i<Workspace::wsv_data.nelem(); ++i )
  {
    if (Workspace::wsv_data[i].Group() == get_wsv_group_id("Agenda"))
    {
      os << "<li><a href=\"/agendas/" << Workspace::wsv_data[i].Name() << "\">"
      << Workspace::wsv_data[i].Name() << "</a>" << endl;
      hitcount2++;
      
      if (hitcount2 == hitcount/2)
        os << "</ul></td><td><ul>" << endl;
    }
  }
  
  os << "</ul></td></tr>" << endl;
}

void ds_list_groups (ostream &os)
{
  extern const ArrayOfString wsv_group_names;
  Index i;
  
  os << "<tr><td colspan=\"2\"><a name=\"groups\" />"
  << "<h2>Workspace Groups</h2></td></tr>" << endl
  << "<tr><td><ul>" << endl;
  
  for ( i=0; i<wsv_group_names.nelem(); ++i )
  {
    os << "<li><a href=\"/groups/" << wsv_group_names[i] << "\">"
    << wsv_group_names[i] << "</a>" << endl;
    
    if (i == wsv_group_names.nelem()/2)
      os << "</ul></td><td><ul>" << endl;
  }
  
  os << "</ul></td></tr>" << endl;
}

void ds_list_methods (ostream &os)
{
  extern const Array<MdRecord> md_data_raw;
  Index i;
  
  os << "<tr><td colspan=\"2\"><a name=\"methods\" />"
  << "<h2>Workspace Methods</h2></td></tr>" << endl
  << "<tr><td><ul>" << endl;
  
  for ( i=0; i<md_data_raw.nelem(); ++i )
  {
    os << "<li><a href=\"/methods/" << md_data_raw[i].Name() << "\">"
    << md_data_raw[i].Name() << "</a>" << endl;
    
    if (i == md_data_raw.nelem()/2)
      os << "</ul></td><td><ul>" << endl;
  }
  
  os << "</ul></td></tr>" << endl;
}

void ds_list_variables (ostream &os)
{
  Index i;
  
  os << "<tr><td colspan=\"2\"><a name=\"variables\" />"
  << "<h2>Workspace Variables</h2></td></tr>" << endl
  << "<tr><td><ul>" << endl;
  
  Index hitcount = 0;
  for ( i=0; i<Workspace::wsv_data.nelem(); ++i )
  {
    if (Workspace::wsv_data[i].Group() != get_wsv_group_id("Agenda"))
      hitcount++;
  }

  Index hitcount2 = 0;
  for ( i=0; i<Workspace::wsv_data.nelem(); ++i )
  {
    if (Workspace::wsv_data[i].Group() != get_wsv_group_id("Agenda"))
    {
      os << "<li>" << ds_insert_wsv_link(Workspace::wsv_data[i].Name()) << endl;
      hitcount2++;
      
      if (hitcount2 == hitcount/2)
        os << "</ul></td><td><ul>" << endl;
    }
  }
  
  os << "</ul></td></tr>" << endl;
}

void ds_doc_method (ostream &os, const string& mname)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data_raw;
  extern const map<String, Index> MdRawMap;
  extern const ArrayOfString wsv_group_names;
  
  // Let's first assume it is a method that the user wants to have
  // described.
  
  // Find method id:
  map<String, Index>::const_iterator it =
  MdRawMap.find(mname);
  if ( it != MdRawMap.end() )
  {
    // If we are here, then the given name matches a method.
    const MdRecord& mdr = md_data_raw[it->second];
    String indent = "";
    
    os << "<h3>Description</h3>" << endl;
    
    char ch = 0;
    for (String::const_iterator sit = mdr.Description().begin();
         sit != mdr.Description().end(); sit++)
    {
      os << *sit;
      if (ch == '\n' && *sit == '\n') os << "<p>";
      ch = *sit;
    }

    os << endl;
    
    bool is_first_author = true;
    for (Index i = 0; i < mdr.Authors().nelem(); i++)
    {
      if (is_first_author)
      {
        os << "<p><b>Authors: </b>";
        is_first_author = false;
      }
      else
        os << ", ";
      
      os << mdr.Authors()[i];
    }
    os << "\n";
    
    // Print the method's synopsis
    while (indent.length() < mdr.Name().length() + 2) indent += ' ';
    
    os << "<h3>Synopsis</h3>" << endl;
    
    ostringstream buf;
    ostringstream param;
    const size_t linelen = 2048;
    
    buf << "<p><table><tr><td>" << mdr.Name() << "(&nbsp;</td><td>";
    bool first = true;
    for ( Index i=0; i<mdr.Out().nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      param << ds_insert_wsv_link(Workspace::wsv_data[mdr.Out()[i]].Name());
      
      limit_line_length( os, buf, param, indent, linelen );
    }
    
    for ( Index i=0; i<mdr.GOutType().nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      if (mdr.GOut()[i].length())
        param << mdr.GOut()[i];
      else
        param << "gout" << i;
      
      limit_line_length( os, buf, param, indent, linelen );
    }
    
    const ArrayOfIndex &inonly = mdr.InOnly();
    for ( Index i=0; i<inonly.nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      param << ds_insert_wsv_link(Workspace::wsv_data[inonly[i]].Name());
      
      limit_line_length( os, buf, param, indent, linelen );
    }
    
    for ( Index i=0; i<mdr.GInType().nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      if (mdr.GIn()[i].length())
      {
        param << mdr.GIn()[i];
      }
      else
      {
        param << "gin" << i;
      }
      
      limit_line_length( os, buf, param, indent, linelen );
    }
    if (buf.str().length()) os << buf.str();
    
    os << " )</td></tr></table>" << endl;
    
    os << "<h3>Variables</h3>" << endl;
    
    // Out:
    indent = "";
    String desc;
    os << "<table>" << endl;
    for ( Index i=0; i<mdr.Out().nelem(); ++i )
    {
      buf.str("");
      buf << "<tr>";
      buf <<    "<td>OUT</td>";
      
      {
        const String& vname = Workspace::wsv_data[mdr.Out()[i]].Name();
        buf << "<td align=\"right\">" << ds_insert_wsv_link(vname) << "</td><td>(";
        buf << ds_insert_group_link(wsv_group_names[Workspace::wsv_data[mdr.Out()[i]].Group()]);
        buf << ")</td><td>";
      }
      
      get_short_wsv_description(desc, Workspace::wsv_data[mdr.Out()[i]].Description());
      
      if (buf.str().length() + desc.length() > linelen)
      {
        format_paragraph (desc, indent, linelen);
        buf << endl << indent << desc;
      }
      else
      {
        buf << desc;
      }
      
      os << buf.str() << "</td></tr>" << endl;
    }
    
    size_t lastlen;
    bool fit;
    for ( Index i=0; i<mdr.GOut().nelem(); ++i )
    {
      buf.str("");
      buf << "<tr>";
      buf <<    "<td>GOUT</td><td align=\"right\">" << mdr.GOut()[i] << "</td><td>(";
      if (mdr.GOutType()[i] == get_wsv_group_id("Any")
          && mdr.GOutSpecType()[i].nelem())
      {
        bool firstarg = true;
        for (Index j = 0; j < mdr.GOutSpecType()[i].nelem(); j++)
        {
          if (!firstarg) buf << ", "; else firstarg = false;
          buf << ds_insert_group_link(wsv_group_names[mdr.GOutSpecType()[i][j]]);
        }
      }
      else
      {
        buf << ds_insert_group_link(wsv_group_names[mdr.GOutType()[i]]);
      }
      
      buf << ")</td><td>";
      desc = buf.str();
      lastlen = desc.length();
      fit = format_paragraph (desc, indent, linelen);
      buf.str("");
      os << desc;
      
      desc = mdr.GOutDescription()[i];
      if (!fit)
      {
        format_paragraph (desc, indent, linelen);
        buf << endl << indent << desc;
      }
      else if (lastlen + desc.length() > linelen)
      {
        format_paragraph (desc, indent, linelen, lastlen);
        buf << endl << desc;
      }
      else
      {
        buf << desc;
      }
      
      os << buf.str() << "</td></tr>" << endl;
    }
    
    for ( Index i=0; i<mdr.In().nelem(); ++i )
    {
      buf.str("");
      buf << "<tr>";
      buf <<    "<td>IN</td>";
      
      const String& vname = Workspace::wsv_data[mdr.In()[i]].Name();
      buf << "<td align=\"right\">" << ds_insert_wsv_link(vname);
      buf << "</td><td>(";
      buf << ds_insert_group_link(wsv_group_names[Workspace::wsv_data[mdr.In()[i]].Group()]);
      buf << ")</td><td>";
      
      get_short_wsv_description(
                                desc, Workspace::wsv_data[mdr.In()[i]].Description());
      
      if (buf.str().length() + desc.length() > linelen)
      {
        format_paragraph (desc, indent, linelen, indent.length());
        buf << endl << indent << desc;
      }
      else
      {
        buf << desc;
      }
      
      os << buf.str() << "</td></tr>" << endl;
    }
    
    for ( Index i=0; i<mdr.GIn().nelem(); ++i )
    {
      buf.str("");
      buf << "<tr>";
      buf <<    "<td>GIN</td><td align=\"right\">" << mdr.GIn()[i] << "</td><td>(";
      if (mdr.GInType()[i] == get_wsv_group_id("Any")
          && mdr.GInSpecType()[i].nelem())
      {
        bool firstarg = true;
        for (Index j = 0; j < mdr.GInSpecType()[i].nelem(); j++)
        {
          if (!firstarg) buf << ", "; else firstarg = false;
          buf << ds_insert_group_link(wsv_group_names[mdr.GInSpecType()[i][j]]);
        }
      }
      else
      {
        buf << ds_insert_group_link(wsv_group_names[mdr.GInType()[i]]);
      }
      
      if (mdr.GInDefault()[i] != NODEF)
      {
        buf << ", Default: ";
        if (mdr.GInType()[i] == get_wsv_group_id ("String"))
        {
          buf << "\"" << mdr.GInDefault()[i] << "\"";
        }
        else
        {
          buf << mdr.GInDefault()[i];
        }
        
      }
      
      buf << ")</td><td>";
      desc = buf.str();
      lastlen = desc.length();
      fit = format_paragraph (desc, indent, linelen);
      buf.str("");
      os << desc;
      
      desc = mdr.GInDescription()[i];
      if (!fit)
      {
        format_paragraph (desc, indent, linelen);
        buf << indent << desc;
      }
      else if (lastlen + desc.length() > linelen)
      {
        format_paragraph (desc, indent, linelen, indent.length());
        buf << indent << desc;
      }
      else
      {
        buf << desc;
      }
      
      os << buf.str() << "</td></tr>" << endl;
    }
    os << "</table>" << endl;
  }
}

void ds_doc_variable_methods(ostream& os, const string& vname)
{
  // Check if the user gave the name of a specific variable.
  map<String, Index>::const_iterator mi =
  Workspace::WsvMap.find(vname);
  extern const Array<MdRecord>  md_data_raw;
  if ( mi != Workspace::WsvMap.end() )
  {
    // If we are here, then the given name matches a variable.
    Index wsv_key = mi->second;
    
    // List generic methods:
    Index hitcount = 0;
    os << "<h3>Generic and supergeneric methods that can generate " << vname << "</h3>" << endl;
    os << "<ul>" << endl;
    for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];
      
      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if ( count( mdd.GOutType().begin(),
                 mdd.GOutType().end(),
                 Workspace::wsv_data[wsv_key].Group() ) )
      {
        os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      }
      else if  ( count( mdd.GOutType().begin(),
                       mdd.GOutType().end(),
                       get_wsv_group_id("Any") ) )
      {
        for (Index j = 0; j < mdd.GOutType().nelem(); j++)
        {
          if (mdd.GOutType()[j] == get_wsv_group_id("Any"))
          {
            if (mdd.GOutSpecType()[j].nelem())
            {
              if (count( mdd.GOutSpecType()[j].begin(),
                        mdd.GOutSpecType()[j].end(),
                        Workspace::wsv_data[wsv_key].Group() ) )
              {
                os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            }
            else
            {
              os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if ( 0==hitcount ) os << "<li>none" << endl;
    
    os << endl << "</ul>" << endl;
    
    // List specific methods:
    hitcount = 0;
    os 
    << "<h3>Specific methods that can generate " << vname << "</h3>" << endl
    << "<ul>" << endl;
    for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];
      
      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if ( count( mdd.Out().begin(),
                 mdd.Out().end(),
                 wsv_key ) ) 
      {
        os << "<li>" << ds_insert_wsm_link(mdd.Name()) << "\n";
        ++hitcount;
      }
    }
    if ( 0==hitcount ) os << "<li>none\n";
    
    os << endl << "</ul>" << endl;
    
    // List generic methods:
    hitcount = 0;
    os << "<h3>Generic and supergeneric methods that can use " << vname << "</h3>" << endl;
    os << "<ul>" << endl;
    for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];
      
      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if ( count( mdd.GInType().begin(),
                 mdd.GInType().end(),
                 Workspace::wsv_data[wsv_key].Group() ) )
      {
        os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      }
      else if  ( count( mdd.GInType().begin(),
                       mdd.GInType().end(),
                       get_wsv_group_id("Any") ) )
      {
        for (Index j = 0; j < mdd.GInType().nelem(); j++)
        {
          if (mdd.GInType()[j] == get_wsv_group_id("Any"))
          {
            if (mdd.GInSpecType()[j].nelem())
            {
              if (count( mdd.GInSpecType()[j].begin(),
                        mdd.GInSpecType()[j].end(),
                        Workspace::wsv_data[wsv_key].Group() ) )
              {
                os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            }
            else
            {
              os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if ( 0==hitcount ) os << "<li>none" << endl;
    
    os << endl << "</ul>" << endl;
    
    // List specific methods:
    hitcount = 0;
    os 
    << "<h3>Specific methods that require " << vname << "</h3>" << endl
    << "<ul>" << endl;
    for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];
      
      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if ( count( mdd.In().begin(),
                 mdd.In().end(),
                 wsv_key ) ) 
      {
        os << "<li>" << ds_insert_wsm_link(mdd.Name()) << "\n";
        ++hitcount;
      }
    }
    if ( 0==hitcount ) os << "<li>none\n";
    
    os << endl << "</ul>" << endl;
    
    // List agendas with this variable as output:
    extern Array<AgRecord> agenda_data;
    hitcount = 0;
    os 
    << "<h3>Agendas that can generate " << vname << "</h3>" << endl
    << "<ul>" << endl;
    for ( Index i=0; i<agenda_data.nelem(); ++i )
    {
      // Get handle on method record:
      const AgRecord& ar = agenda_data[i];
      
      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if ( count( ar.Out().begin(),
                 ar.Out().end(),
                 wsv_key ) ) 
      {
        os << "<li>" << ds_insert_agenda_link(ar.Name()) << "\n";
        ++hitcount;
      }
    }
    if ( 0==hitcount ) os << "<li>none\n";
    
    os << endl << "</ul>" << endl;
    
    // List agendas with this variable as input:
    hitcount = 0;
    os << "<h3>Agendas that require " << vname << "</h3>" << endl
    << "<ul>" << endl;
    for ( Index i=0; i<agenda_data.nelem(); ++i )
    {
      // Get handle on method record:
      const AgRecord& ar = agenda_data[i];
      
      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if ( count( ar.In().begin(),
                 ar.In().end(),
                 wsv_key ) ) 
      {
        os << "<li>" << ds_insert_agenda_link(ar.Name()) << "\n";
        ++hitcount;
      }
    }
    
    if ( 0==hitcount ) os << "<li>none\n";
    
    os << endl << "</ul>" << endl;
  }
}

void ds_doc_variable (ostream &os, const string& vname)
{
  extern const ArrayOfString wsv_group_names;
  
  // Find wsv id:
  map<String, Index>::const_iterator it = Workspace::WsvMap.find(vname);
  if ( it != Workspace::WsvMap.end() )
  {
    // If we are here, then the given name matches a workspace
    // variable.
    char ch = 0;
    for (String::const_iterator sit = Workspace::wsv_data[it->second].Description().begin();
         sit != Workspace::wsv_data[it->second].Description().end(); sit++)
    {
      os << *sit;
      if (ch == '\n' && *sit == '\n') os << "<p>";
      ch = *sit;
    }
    
    os << endl;
    
    os << "<p><b>Group: </b>"
    << ds_insert_group_link(wsv_group_names[Workspace::wsv_data[it->second].Group()]) << endl;
  }

  ds_doc_variable_methods(os, vname);
}

void ds_doc_agenda (ostream &os, const string& aname)
{
  extern const ArrayOfString wsv_group_names;
  
  // Find wsv id:
  map<String, Index>::const_iterator it = Workspace::WsvMap.find(aname);
  extern Array<AgRecord> agenda_data;
  extern const map<String, Index> AgendaMap;
  map<String, Index>::const_iterator ait = AgendaMap.find(aname);
  
  if ( it != Workspace::WsvMap.end() && ait != AgendaMap.end() )
  {
    // If we are here, then the given name matches a workspace
    // variable.
    char ch = 0;
    for (String::const_iterator sit = agenda_data[ait->second].Description().begin();
         sit != agenda_data[ait->second].Description().end(); sit++)
    {
      os << *sit;
      if (ch == '\n' && *sit == '\n') os << "<p>";
      ch = *sit;
    }
    
    os << endl;
    
    os << "<p><b>Group: </b>"
    << ds_insert_group_link(wsv_group_names[Workspace::wsv_data[it->second].Group()]) << endl;
    
    os << "<h3>Variables</h3>" << endl;
    
    // Out:
    if ( ait != AgendaMap.end() )
    {
      // If we are here, then the given name matches a method.
      const AgRecord& agr = agenda_data[ait->second];
      String indent = "";
      String desc;
      ostringstream buf;
      size_t linelen = 80;
      os << "<table>" << endl;
      for ( Index i=0; i<agr.Out().nelem(); ++i )
      {
        buf.str("");
        buf << "<tr>";
        buf <<    "<td>OUT</td>";
        
        {
          const String& vname = Workspace::wsv_data[agr.Out()[i]].Name();
          buf << "<td align=\"right\">" << ds_insert_wsv_link(vname) << "</td><td>(";
          buf << ds_insert_group_link(wsv_group_names[Workspace::wsv_data[agr.Out()[i]].Group()]);
          buf << ")</td><td>";
        }
        
        get_short_wsv_description(desc, Workspace::wsv_data[agr.Out()[i]].Description());
        
        if (buf.str().length() + desc.length() > linelen)
        {
          format_paragraph (desc, indent, linelen);
          buf << endl << indent << desc;
        }
        else
        {
          buf << desc;
        }
        
        os << buf.str() << "</td></tr>" << endl;
        
      }
      
      for ( Index i=0; i<agr.In().nelem(); ++i )
      {
        buf.str("");
        buf << "<tr>";
        buf <<    "<td>IN</td>";
        
        const String& vname = Workspace::wsv_data[agr.In()[i]].Name();
        buf << "<td align=\"right\">" << ds_insert_wsv_link(vname);
        buf << "</td><td>(";
        buf << ds_insert_group_link(wsv_group_names[Workspace::wsv_data[agr.In()[i]].Group()]);
        buf << ")</td><td>";
        
        get_short_wsv_description(
                                  desc, Workspace::wsv_data[agr.In()[i]].Description());
        
        if (buf.str().length() + desc.length() > linelen)
        {
          format_paragraph (desc, indent, linelen, indent.length());
          buf << endl << indent << desc;
        }
        else
        {
          buf << desc;
        }
        
        os << buf.str() << "</td></tr>" << endl;
      }
      
      os << "</table>" << endl;
    }
  }

  ds_doc_variable_methods(os, aname);
}

void ds_doc_group (ostream &os, const string& gname)
{
  // Check if the user gave the name of a specific variable.
  Index gid = get_wsv_group_id (gname);
  extern const Array<MdRecord>  md_data_raw;
  if ( gid != -1 )
  {
    // If we are here, then the given name matches a group.

    // List generic methods:
    Index hitcount = 0;
    os << "<h3>Generic and supergeneric methods that can generate " << gname << "</h3>" << endl;
    os << "<ul>" << endl;
    for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];
      
      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if ( count( mdd.GOutType().begin(), mdd.GOutType().end(), gid ) )
      {
        os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      }
      else if  ( count( mdd.GOutType().begin(),
                       mdd.GOutType().end(),
                       get_wsv_group_id("Any") ) )
      {
        for (Index j = 0; j < mdd.GOutType().nelem(); j++)
        {
          if (mdd.GOutType()[j] == get_wsv_group_id("Any"))
          {
            if (mdd.GOutSpecType()[j].nelem())
            {
              if (count( mdd.GOutSpecType()[j].begin(),
                        mdd.GOutSpecType()[j].end(),
                        gid ) )
              {
                os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            }
            else
            {
              os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if ( 0==hitcount ) os << "<li>none" << endl;
    
    os << endl << "</ul>" << endl;

    hitcount = 0;
    os << "<h3>Generic and supergeneric methods that can use " << gname << "</h3>" << endl;
    os << "<ul>" << endl;
    for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];
      
      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if ( count( mdd.GInType().begin(), mdd.GInType().end(), gid ) )
      {
        os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      }
      else if  ( count( mdd.GInType().begin(),
                       mdd.GInType().end(),
                       get_wsv_group_id("Any") ) )
      {
        for (Index j = 0; j < mdd.GInType().nelem(); j++)
        {
          if (mdd.GInType()[j] == get_wsv_group_id("Any"))
          {
            if (mdd.GInSpecType()[j].nelem())
            {
              if (count( mdd.GInSpecType()[j].begin(),
                        mdd.GInSpecType()[j].end(),
                        gid ) )
              {
                os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            }
            else
            {
              os << "<li>" << ds_insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if ( 0==hitcount ) os << "<li>none" << endl;
    
    os << endl << "</ul>" << endl;

    if (gname != "Any")
    {
      // List specific methods:
      hitcount = 0;
      os << "<h3>Specific methods that can generate " << gname << "</h3>" << endl;
      os << "<ul>" << endl;
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
      {
        // Get handle on method record:
        const MdRecord& mdd = md_data_raw[i];
        
        bool first = true;
        for ( Index j=0; j<mdd.Out().nelem(); j++)
        {
          // This if statement checks whether the type of this output variable
          // matches this group.
          if (Workspace::wsv_data[mdd.Out()[j]].Group() == gid)
          {
            if (first)
            {
              first = false;
              os << "<li>" << ds_insert_wsm_link(mdd.Name()) << " (";
            }
            else
              os << ", ";
            os << ds_insert_wsv_link(Workspace::wsv_data[mdd.Out()[j]].Name());
            
            ++hitcount;
          }
        }
        if (!first) os << ")" << endl;
      }
      if ( 0==hitcount ) os << "<li>none" << endl;
      
      os << endl << "</ul>" << endl;
      
      hitcount = 0;
      os << "<h3>Specific methods that require variables of group " << gname << "</h3>" << endl;
      os << "<ul>" << endl;
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
      {
        // Get handle on method record:
        const MdRecord& mdd = md_data_raw[i];
        
        bool first = true;
        for ( Index j=0; j<mdd.In().nelem(); j++)
        {
          // This if statement checks whether the type of this output variable
          // matches this group.
          if (Workspace::wsv_data[mdd.In()[j]].Group() == gid)
          {
            if (first)
            {
              first = false;
              os << "<li>" << ds_insert_wsm_link(mdd.Name()) << " (";
            }
            else
              os << ", ";
            os << ds_insert_wsv_link(Workspace::wsv_data[mdd.In()[j]].Name());
            
            ++hitcount;
          }
        }
        if (!first) os << ")" << endl;
      }
      if ( 0==hitcount ) os << "<li>none" << endl;
      
      os << endl << "</ul>" << endl;
      
      Index i;
      
      os << "<h3>Workspace Variables of group " << gname << "</h3>" << endl
      << "<ul>" << endl;
      
      hitcount = 0;
      for ( i=0; i<Workspace::wsv_data.nelem(); ++i )
      {
        if (Workspace::wsv_data[i].Group() == get_wsv_group_id(gname))
        {
          os << "<li>" << ds_insert_wsv_link(Workspace::wsv_data[i].Name()) << endl;
          hitcount++;
        }
      }
      if ( 0==hitcount ) os << "<li>none" << endl;
      
      os << "</ul></td></tr>" << endl;
    }
  }
}

void ds_insert_breadcrumb_token (ostream& os, const vector<string>& tokens,
                                 size_t token_id)
{
  bool link = (token_id < tokens.size()-1);
  
  if (link)
  {
    os << "<a href=\"";
    for (size_t t = 0; t <= token_id; t++)
    {
      if (t != 1) os << "/";
      os << tokens[t];
    }
    os << "\">";
  }

  if (!token_id) os << "Home";
  else if (tokens[token_id] == "methods") os << "Methods";
  else if (tokens[token_id] == "variables") os << "Variables";
  else if (tokens[token_id] == "agendas") os << "Agendas";
  else if (tokens[token_id] == "groups") os << "Groups";
  else os << tokens[token_id];
  
  if (link) os << "</a>";
}

void ds_insert_breadcrumbs (ostream& os, const vector<string>& tokens)
{
  os << "<span id=\"breadcrumbs\">";
  for (size_t t = 0; t < tokens.size(); t++)
  {
    if (t) os << "&nbsp;>>&nbsp;";
    ds_insert_breadcrumb_token (os, tokens, t);
  }
  os << "</span>" << endl;
}

void ds_insert_title (ostream& os, const string& title = "")
{
  os << "<h1>"DOCSERVER_NAME;
  if (title.length())
    os << " &mdash; " << title;
  os << "</h1>" << endl;
}

void ds_insert_index (ostream& os, const vector<string>& tokens)
{
  void(*index_method)(ostream&);
  string title;
  
  if (tokens.size() < 2)
  {
    ds_insert_breadcrumbs (os, tokens);
    ds_insert_title (os, "Index");
    
    os
    << "<p>Jump to: <a href=\"#methods\">Methods</a>&nbsp;-&nbsp;"
    << "<a href=\"#variables\">Variables</a>&nbsp;-&nbsp;"
    << "<a href=\"#agendas\">Agendas</a>&nbsp;-&nbsp;"
    << "<a href=\"#groups\">Groups</a>"
    << endl;
    
    os << "<center><table width=\"90%\">" << endl;
    ds_list_methods (os);
    ds_list_variables (os);
    ds_list_agendas (os);
    ds_list_groups (os);
    os << "</table></center>" << endl;
    return;
  }
  else if (tokens[1] == "methods")
  {
    title = "Workspace Method Index";
    index_method = ds_list_methods;
  }
  else if (tokens[1] == "variables")
  {
    title = "Workspace Variable Index";
    index_method = ds_list_variables;
  }
  else if (tokens[1] == "groups")
  {
    title = "Workspace Group Index";
    index_method = ds_list_groups;
  }
  else if (tokens[1] == "agendas")
  {
    title = "Agenda Index";
    index_method = ds_list_agendas;
  }
  else return;

  ds_insert_breadcrumbs (os, tokens);
  ds_insert_title (os, title);
  os << "<center><table width=\"90%\">" << endl;
  (*index_method)(os);
  os << "</table></center>" << endl;
}

void ds_insert_doc (ostream& os, const vector<string>& tokens)
{
  void(*doc_method)(ostream&, const string&);
  string title;
  
  if (tokens[1] == "methods")
  {
    title = "Workspace Method " + tokens[2];
    doc_method = ds_doc_method;
  }
  else if (tokens[1] == "variables")
  {
    title = "Workspace Variable " + tokens[2];
    doc_method = ds_doc_variable;
  }
  else if (tokens[1] == "groups")
  {
    title = "Workspace Group " + tokens[2];
    doc_method = ds_doc_group;
  }
  else if (tokens[1] == "agendas")
   {
   title = "Agenda " + tokens[2];
   doc_method = ds_doc_agenda;
   }
  else return;
  
  ds_insert_breadcrumbs (os, tokens);
  ds_insert_title (os);
  os << "<h2>" << title << "</h2>" << endl;
  
  (*doc_method)(os, tokens[2]);
}

void ds_stylesheet (ostream& os)
{
  os
  << "body { font-family: monospace }" << endl
  << "a:link { color: #3465a4; text-decoration: underline; }" << endl
  << "a:visited { color: #729fcf; text-decoration: underline; }" << endl
  << "a:active { color: #ce5c00; text-decoration: underline; background-color: #eeeeec}" << endl
  << "a:hover { color: #f57900; text-decoration: underline; }" << endl
  
  << "table {" << endl
  << "border: 0px;" << endl
  << "}" << endl
  
  << "table td {" << endl
  << "vertical-align: top;" << endl
  << "}" << endl
  
  << "span#breadcrumbs {" << endl
  << "font-size: small;" << endl
  << "}" << endl
  
  << endl;
}

static int
ahc_echo (void *cls _U_,
          struct MHD_Connection *connection,
          const char *url,
          const char *method,
          const char *version _U_,
          const char *upload_data _U_, size_t *upload_data_size _U_, void **ptr)
{
  static int aptr;
  const char *val;
  string surl(url);
  struct MHD_Response *response;
  int ret;
  
  if (0 != strcmp (method, "GET"))
  {
    out0 << "Docserver error: Unexpected method " << method << "\n";
    return MHD_NO;              /* unexpected method */
  }
  if (&aptr != *ptr)
  {
    /* do never respond on first call */
    *ptr = &aptr;
    return MHD_YES;
  }
  *ptr = NULL;                  /* reset when done */
  val = MHD_lookup_connection_value (connection, MHD_GET_ARGUMENT_KIND, "q");
  
  ostringstream hout;
  
  vector<string> tokens;
  split (surl, '/', tokens);
  
  string content_type = "text/html; charset=utf-8";
  if (tokens.size() == 2 && tokens[1] == "styles.css")
  {
    ds_stylesheet(hout);
    content_type = "text/css; charset=utf-8";
  }
  else
  {
    ds_begin_page (hout);
    
    switch (tokens.size())
    {
      case 0:
      case 1:
      case 2:
        ds_insert_index (hout, tokens);
        break;
      case 3:
        ds_insert_doc(hout, tokens);
        break;
      default:
        ds_insert_title (hout);
    }
    
    ds_end_page (hout);
  }
  
  response = MHD_create_response_from_data (hout.str().length(),
                                            (void *)hout.str().c_str(),
                                            MHD_NO, MHD_YES);
  
  if (response == NULL) {
    out0 << "Docserver error: response = 0\n";
    return MHD_NO;
  }

  MHD_add_response_header (response, "Content-type", content_type.c_str());
  ret = MHD_queue_response (connection, MHD_HTTP_OK, response);
  MHD_destroy_response (response);
  
  return ret;
}

int docserver_start(Index port, bool daemon)
{
  struct MHD_Daemon *d;
  
  if (port == -1) port=9000;
  d = MHD_start_daemon (MHD_USE_THREAD_PER_CONNECTION | MHD_USE_DEBUG,
                        port,
                        NULL, NULL, &ahc_echo, NULL, MHD_OPTION_END);
  if (d == NULL)
  {
    out0 << "Error: Cannot start server. Maybe port " << port << " is already in use?\n"; 
    return 1;
  }
  else
  {
    if (daemon)
      out0 << "ARTS docserver listening at http://localhost:" << port << "\n";
    else
      out0 << "\n"
      << "===========================================================\n"
      << "Now point your web browser to http://localhost:" << port << "\n"
      << "===========================================================\n\n"
      << "Press enter to exit.\n";
  }
  
  
  if (daemon)
  {
    pause();
  }
  else
  {
    (void) getc (stdin);
    out0 << "Stopping docserver.\n";
    MHD_stop_daemon (d);
    out0 << "Goodbye.\n";
  }
  return 0;
}

