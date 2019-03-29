<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"     xmlns:xsl="http://www.w3.org/1999/XSL/Transform"    xmlns:mml="http://morphml.org/morphml/schema"    xmlns:meta="http://morphml.org/metadata/schema"    xmlns:nml="http://morphml.org/neuroml/schema"    xmlns:cml="http://morphml.org/channelml/schema"    xmlns:bio="http://morphml.org/biophysics/schema"    xmlns:net="http://morphml.org/networkml/schema"    exclude-result-prefixes="mml meta nml net cml bio">

<!--

    This file is used to convert NeuroML files (morphology and/or network structure)
    to X3D files, for visualisation of 3D structure in any browser with an X3D plugin or 
    X3D standalone application
    
    Funding for this work has been received from the Medical Research Council and the 
    Wellcome Trust. This file was initially developed as part of the neuroConstruct project
    
    Author: Padraig Gleeson
    Copyright 2009 University College London
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
-->

<xsl:variable name="showAxes">0</xsl:variable>
<xsl:variable name="defaultCellRadius">5</xsl:variable>

<xsl:output method="xml" indent="yes" />

<xsl:template  match="/">
<X3D profile="Immersive.." version="2.0">
    
    <Scene>
        <Background   skyColor="0.6 0.7 0.9"/>
        <Viewpoint description="Down z axis, 500 microns away" position="0 0 500"/> 
        <Viewpoint description="Down z axis, 200 microns away" position="0 0 200"/> 
        <Viewpoint description="Down z axis, 2mm away" position="0 0 2000"/> 
        
    <xsl:if test="$showAxes = 1">

        <xsl:call-template name="showAxes"/>
    
    </xsl:if> 
    
    <xsl:apply-templates select="nml:neuroml"/>
    
    <xsl:apply-templates select="mml:morphml"/>
    
    <xsl:apply-templates select="net:networkml"/>
    
    </Scene>
    
</X3D>

</xsl:template>

<xsl:template match="mml:morphml | nml:neuroml">
   
    <xsl:for-each select="mml:cells/mml:cell | nml:cells/nml:cell">
        
        <xsl:for-each select="mml:segments/mml:segment">
            
            <xsl:variable name="proximal">
                <xsl:choose>
                    <xsl:when test="count(mml:proximal) &gt; 0"> 
                        <xsl:value-of select="mml:proximal/@x"/><xsl:text>  </xsl:text><xsl:value-of select="mml:proximal/@y"/><xsl:text>  </xsl:text><xsl:value-of select="mml:proximal/@z"/>
                        </xsl:when>
                    <xsl:otherwise><xsl:variable name="parent"><xsl:value-of select="@parent"/></xsl:variable>
                    <xsl:for-each select="../mml:segment[@id = $parent]"><xsl:value-of select="mml:distal/@x"/><xsl:text>  </xsl:text><xsl:value-of select="mml:distal/@y"/><xsl:text>  </xsl:text><xsl:value-of select="mml:distal/@z"/></xsl:for-each></xsl:otherwise>
                </xsl:choose>
            </xsl:variable>
                
            <xsl:variable name="distal">
                <xsl:value-of select="mml:distal/@x"/><xsl:text>  </xsl:text><xsl:value-of select="mml:distal/@y"/><xsl:text>  </xsl:text><xsl:value-of select="mml:distal/@z"/>
            </xsl:variable>


            <Transform>
                <Shape>
                    <Appearance>
                        <Material/>
                    </Appearance>
                    <LineSet vertexCount="2">
                        <xsl:element name="Coordinate">
                            <xsl:attribute name="point"><xsl:value-of select="$proximal"/>, <xsl:value-of select="$distal"/></xsl:attribute>
                        </xsl:element>
                        <Color color="0 0 0, 0 0 0"/>
                    </LineSet>
                </Shape>
            </Transform>
            <!--
            <xsl:element name="Transform">
                <xsl:attribute name="translation"><xsl:value-of select="$proximal"/></xsl:attribute>
                <Shape>
                    <Appearance>
                      <Material diffuseColor="0 1 0"/>
                    </Appearance>
                    <xsl:element name="Sphere">
                        <xsl:attribute name="radius"><xsl:value-of select="number(mml:proximal/@diameter)*0.5"/></xsl:attribute>
                    </xsl:element>
                </Shape>
            </xsl:element>
            
            <xsl:element name="Transform">
                <xsl:attribute name="translation"><xsl:value-of select="$distal"/></xsl:attribute>
                <Shape>
                    <Appearance>
                      <Material diffuseColor="0 1 0"/>
                    </Appearance>
                    <xsl:element name="Sphere">
                        <xsl:attribute name="radius"><xsl:value-of select="number(mml:distal/@diameter)*0.5"/></xsl:attribute>
                    </xsl:element>
                </Shape>
            </xsl:element>-->
            
            
        </xsl:for-each>
        
        
    </xsl:for-each>
    
    
</xsl:template>

<xsl:template match="net:networkml">
        
    <xsl:for-each select="net:populations/net:population">
        <xsl:variable name="cell_type" select="@cell_type"/>
        <xsl:for-each select="net:instances/net:instance">
            
            <xsl:variable name="location"><xsl:value-of select="net:location/@x"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@y"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@z"/></xsl:variable>
            <xsl:element name="Transform">
                <xsl:attribute name="translation"><xsl:value-of select="$location"/></xsl:attribute>
            <xsl:choose>
              <xsl:when test="$cell_type='mitral'">
        <Shape>
            <Appearance>
              <Material diffuseColor="1 0 0"/>
            </Appearance>
            <xsl:element name="Sphere">
                <xsl:attribute name="radius"><xsl:value-of select="0.00003"/></xsl:attribute>
            </xsl:element>
        </Shape>
              </xsl:when>
              <xsl:otherwise>
        <Shape>
            <Appearance>
              <Material diffuseColor="0 1 0"/>
            </Appearance>
            <xsl:element name="Sphere">
                <xsl:attribute name="radius"><xsl:value-of select="0.000025"/></xsl:attribute>
            </xsl:element>
        </Shape>
              </xsl:otherwise>
            </xsl:choose>
            </xsl:element>
        </xsl:for-each>
    </xsl:for-each>
        
    <xsl:for-each select="net:projections/net:projection">
        <xsl:variable name="src"><xsl:value-of select="net:source"/><xsl:value-of select="@source"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
        <xsl:variable name="tgt"><xsl:value-of select="net:target"/><xsl:value-of select="@target"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
        
        <xsl:comment>Projection <xsl:value-of select="@name"/> between <xsl:value-of select="$src"/> and <xsl:value-of select="$tgt"/></xsl:comment>
        
        <xsl:for-each select="net:connections/net:connection">
            <xsl:variable name="preCellId"><xsl:value-of select="net:pre/@cell_id"/><xsl:value-of select="@pre_cell_id"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
            <xsl:variable name="postCellId"><xsl:value-of select="net:post/@cell_id"/><xsl:value-of select="@post_cell_id"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
            
            <xsl:variable name="preLocation">
                <xsl:for-each select="../../../../net:populations/net:population[@name = $src]/net:instances/net:instance[@id = $preCellId]">
                    <xsl:value-of select="net:location/@x"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@y"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@z"/>
                </xsl:for-each>
            </xsl:variable>
            
            <xsl:variable name="postLocation">
                <xsl:for-each select="../../../../net:populations/net:population[@name = $tgt]/net:instances/net:instance[@id = $postCellId]">
                    <xsl:value-of select="net:location/@x"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@y"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@z"/>
                </xsl:for-each>
            </xsl:variable>
            
            <Transform>
                <Shape>
                    <Appearance>
                        <Material/>
                    </Appearance>
                    <LineSet vertexCount="2">
                        <xsl:element name="Coordinate">
                            <xsl:attribute name="point"><xsl:value-of select="$preLocation"/>, <xsl:value-of select="$postLocation"/></xsl:attribute>
                        </xsl:element>
                        <Color color="0 1 0, 1 0 0"/><!-- Green to red-->
                    </LineSet>
                </Shape>
            </Transform>
            
            
            
        </xsl:for-each> 
    </xsl:for-each>
    

</xsl:template>

    

<xsl:template name="showAxes">

        <!-- X axis-->
    <Transform rotation="0 0 1 -1.570795">
        <Shape>
            <Appearance><Material diffuseColor="0 1 0"/></Appearance>
            <Cylinder height="200" radius="0.5"/>
        </Shape>
        <Transform translation="0 105 0">
            <Shape>
                <Appearance><Material diffuseColor="0 1 0"/></Appearance>
                <Cone height="10" bottomRadius= "1"/>
            </Shape>
        </Transform>
    </Transform>
    
    <!-- Y axis-->
    <Transform>
        <Shape>
            <Appearance><Material diffuseColor="1 1 0"/></Appearance>
            <Cylinder height="200" radius="0.5"/>
        </Shape>
        <Transform translation="0 105 0">
            <Shape>
                <Appearance><Material diffuseColor="1 1 0"/></Appearance>
                <Cone height="10" bottomRadius= "1"/>
            </Shape>
        </Transform>
    </Transform>
    
    <!-- Z axis-->
    <Transform rotation="1 0 0 1.570795">
        <Shape>
            <Appearance><Material diffuseColor="1 0 0"/></Appearance>
            <Cylinder height="200" radius="0.5"/>
        </Shape>
        <Transform translation="0 105 0">
            <Shape>
                <Appearance><Material diffuseColor="1 0 0"/></Appearance>
                <Cone height="10" bottomRadius= "1"/>
            </Shape>
        </Transform>
    </Transform>
    
</xsl:template>

</xsl:stylesheet>
