<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"     xmlns:xsl="http://www.w3.org/1999/XSL/Transform"    xmlns:mml="http://morphml.org/morphml/schema"    xmlns:meta="http://morphml.org/metadata/schema"    xmlns:nml="http://morphml.org/neuroml/schema"    xmlns:cml="http://morphml.org/channelml/schema"    xmlns:bio="http://morphml.org/biophysics/schema"    xmlns:net="http://morphml.org/networkml/schema"    exclude-result-prefixes="mml meta nml net cml bio">
<!-- NOTE that in the stylesheet below, you must use the namespaces as defined above, not as in the xml document -->
 <!-- trigonometry functions needing for rotation - this module is written by Dimitre Novatchev -->
 <!-- found on the net http://fxsl.sourceforge.net/articles/xslCalculator/The%20FXSL%20Calculator.html -->
 <xsl:import href="trignm.xsl"/>
<!--

    This file is used to convert NeuroML files (morphology and/or network structure)
    to X3D files, for visualisation of 3D structure in any browser with an X3D plugin or 
    X3D standalone application
    
    Funding for this work has been received from the Medical Research Council and the 
    Wellcome Trust. This file was initially developed as part of the neuroConstruct project
    
    Author: Padraig Gleeson
    Copyright 2009 University College London
    
    Author2: Aditya Gilra
    
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
<!-- Set the line width as segment diameter or not?
works only in Octaga; in freewrl, setting linewidthScaleFactor makes blinking, dashed, poorly visible lines.
Now I've set linetype="0", rather than linetype="1" (default) both seem to make 2D flat wide lines, so no loss.
Good thing is that linetype="0" lines do not blink in freewrl.
BUT NOTE linetype="1" works fine on Windows freewrl, though slow/crashing in VBox -->
<!-- true() and false() constructors are needed -->
<xsl:variable name="setLineWidth" select="true()"/>
<xsl:variable name="drawGreyCells" select="true()"/>

<xsl:output method="xml" indent="yes" />

<xsl:template  match="/">
<X3D profile="Immersive" version="3.1" xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.1.xsd '>
    
    <Scene>
        <Background   skyColor="1.0 1.0 1.0"/>
        <!--<Background   skyColor="0.6 0.7 0.9"/>-->
        <!-- Press / key in freewrl to get position and orientation of current view, then negate the first three values of orientation, and use here. -->
        <!-- Press v key to switch between viewpoints. First as needed for schemati figure. -->
        <Viewpoint description="for joints connections"	position="715.9, -558.3, 1306" orientation="-0.732 -0.2904 -0.6163 5.081"/>
        <Viewpoint description="Down z axis, 3mm away" position="0 0 3000"/>
        
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
    <xsl:choose>
      <!-- when populations of cells are present: -->
      <xsl:when test="net:populations">
        <xsl:for-each select="net:populations/net:population">
          <xsl:variable name="cell_type"><xsl:value-of select="@cell_type"/></xsl:variable>
          <xsl:comment>Population <xsl:value-of select="@cell_type"/></xsl:comment>
          <xsl:variable name="colorIndex"><xsl:value-of select="position()"/></xsl:variable>
          <xsl:variable name="popName"><xsl:value-of select="@name"/></xsl:variable>
          <xsl:for-each select="net:instances/net:instance">
              <xsl:variable name="cell_x"><xsl:value-of select="net:location/@x"/></xsl:variable>
              <xsl:variable name="cell_y"><xsl:value-of select="net:location/@y"/></xsl:variable>
              <xsl:variable name="cell_z"><xsl:value-of select="net:location/@z"/></xsl:variable>
              <xsl:variable name="zrotation">
                  <xsl:for-each select="meta:notes">
                      <xsl:if test="contains(text(),'zrotation')"> <!-- does text in the notes contain zrotation? -->
                          <xsl:value-of select="substring-after(.,'zrotation=')"/>
                      </xsl:if>
                  </xsl:for-each>
              </xsl:variable>
              <xsl:variable name="colorRGBA">
                <!--<xsl:call-template name="setColor">
                    <xsl:with-param name="colorIndex" select="$colorIndex"/>
                    <xsl:with-param name="subColorIndex" select="position()"/>                    
                </xsl:call-template>-->
                <xsl:call-template name="setColorByName">
                    <xsl:with-param name="popName" select="$popName"/>
                    <xsl:with-param name="cellId" select="@id"/>                    
                </xsl:call-template>
              </xsl:variable>
              <!--
              <xsl:element name="Transform">
              <xsl:attribute name="translation"><xsl:value-of select="$cell_x"/><xsl:text> </xsl:text><xsl:value-of select="$cell_y"/><xsl:text> </xsl:text><xsl:value-of select="$cell_z"/></xsl:attribute>
                <Shape>
                  <xsl:element name="Text">
                    <xsl:attribute name="string"><xsl:value-of select="@id"/></xsl:attribute>
                    <FontStyle DEF='CenteredFontStyle' justify='"MIDDLE" "MIDDLE"' size="20.0"/>
                  </xsl:element>
                  <Appearance>
                    <Material DEF='DefaultMaterial' diffuseColor='1 1 0'/>
                  </Appearance>
                </Shape>
              </xsl:element>
              -->
              <xsl:if test="$drawGreyCells or not($cell_type='granule') or ((@id=444) or (@id=1160))">
                  <xsl:call-template name="cell_draw">
                    <xsl:with-param name="cell_type" select="$cell_type"/>
                    <xsl:with-param name="x" select="$cell_x"/>
                    <xsl:with-param name="y" select="$cell_y"/>
                    <xsl:with-param name="z" select="$cell_z"/>
                    <!-- normalize-space strips leading and trailing spaces and makes single spaces inbetween -->
                    <xsl:with-param name="zrotation" select="normalize-space($zrotation)"/>
                    <xsl:with-param name="colorRGBA" select="$colorRGBA"/>
                  </xsl:call-template>
              </xsl:if>
          </xsl:for-each>
        </xsl:for-each>
        <xsl:for-each select="net:projections/net:projection">
            <xsl:variable name="src"><xsl:value-of select="net:source"/><xsl:value-of select="@source"/></xsl:variable> <!-- Only one of attr 'source' or sub element <source> should be present-->
            <xsl:variable name="tgt"><xsl:value-of select="net:target"/><xsl:value-of select="@target"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
            <xsl:variable name="projName"><xsl:value-of select="@name"/></xsl:variable> <!-- attr name -->
            <!-- don't draw connections that are just duplicated eg for reciprocal mit-gran connections -->
            <xsl:if test="not( $projName='mitral_granule_extra_exc_joints' or $projName='granule_mitral_inh_joints' 
                    or $projName='mitral_granule_extra_exc_multis' or $projName='granule_mitral_inh_multis'
                    or $projName='mitral_granule_extra_exc_singles' or $projName='granule_mitral_inh_singles'
                    or $projName='PG_mitral' )">
                <xsl:comment>Projection <xsl:value-of select="@name"/> between <xsl:value-of select="$src"/> and <xsl:value-of select="$tgt"/></xsl:comment>
                <xsl:for-each select="net:connections/net:connection">
                    <xsl:variable name="preCellId"><xsl:value-of select="net:pre/@cell_id"/><xsl:value-of select="@pre_cell_id"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
                    <xsl:variable name="postCellId"><xsl:value-of select="net:post/@cell_id"/><xsl:value-of select="@post_cell_id"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
                    <xsl:variable name="preSegmentId"><xsl:value-of select="net:pre/@segment_id"/><xsl:value-of select="@pre_segment_id"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
                    <xsl:variable name="postSegmentId"><xsl:value-of select="net:post/@segment_id"/><xsl:value-of select="@post_segment_id"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
                    <xsl:variable name="pos"><xsl:value-of select="position()-1"/></xsl:variable>
                    <!-- for 1st connection, $pos=-1, but [$pos] will return first element, not error -->
                    <!-- Use preceding-sibling::sibling[1] (used below).
                    It is much faster than indexing by position (commented below)!!!
                    Note that preceding-sibling has reverse axis. So [1] is the just previous element.
                    -->
                    <!--
                    <xsl:variable name="preCellIdOld"><xsl:value-of select="../net:connection[$pos]/net:pre/@cell_id"/><xsl:value-of select="../net:connection[$pos]/@pre_cell_id"/></xsl:variable>
                    <xsl:variable name="postCellIdOld"><xsl:value-of select="../net:connection[$pos]/net:post/@cell_id"/><xsl:value-of select="../net:connection[$pos]/@post_cell_id"/></xsl:variable>
                    <xsl:variable name="preSegmentIdOld"><xsl:value-of select="../net:connection[$pos]/net:pre/@segment_id"/><xsl:value-of select="../net:connection[$pos]/@pre_segment_id"/></xsl:variable>
                    <xsl:variable name="postSegmentIdOld"><xsl:value-of select="../net:connection[$pos]/net:post/@segment_id"/><xsl:value-of select="../net:connection[$pos]/@post_segment_id"/></xsl:variable>
                    -->
                    <xsl:variable name="preCellIdOld"><xsl:value-of select="preceding-sibling::net:connection[1]/net:pre/@cell_id"/><xsl:value-of select="preceding-sibling::net:connection[1]/@pre_cell_id"/></xsl:variable>
                    <xsl:variable name="postCellIdOld"><xsl:value-of select="preceding-sibling::net:connection[1]/net:post/@cell_id"/><xsl:value-of select="preceding-sibling::net:connection[1]/@post_cell_id"/></xsl:variable>
                    <xsl:variable name="preSegmentIdOld"><xsl:value-of select="preceding-sibling::net:connection[1]/net:pre/@segment_id"/><xsl:value-of select="preceding-sibling::net:connection[1]/@pre_segment_id"/></xsl:variable>
                    <xsl:variable name="postSegmentIdOld"><xsl:value-of select="preceding-sibling::net:connection[1]/net:post/@segment_id"/><xsl:value-of select="preceding-sibling::net:connection[1]/@post_segment_id"/></xsl:variable>

                    <!-- draw a connection line only if the connection is not a duplicate of the previous one,
                     or if it is 1st connection -->
                    <xsl:if test="($pos=0) or ($preCellId!=$preCellIdOld) or ($postCellId!=$postCellIdOld) or ($preSegmentId!=$preSegmentIdOld) or ($postSegmentId!=postSegmentIdOld)">
                        <xsl:variable name="preCellType">
                            <xsl:for-each select="/nml:neuroml/net:populations/net:population[@name = $src]">
                                <xsl:value-of select="@cell_type"/>
                            </xsl:for-each>
                        </xsl:variable>
                        <xsl:variable name="postCellType">
                            <xsl:for-each select="/nml:neuroml/net:populations/net:population[@name = $tgt]">
                                <xsl:value-of select="@cell_type"/>
                            </xsl:for-each>
                        </xsl:variable>

                        <xsl:variable name="prePoint">
                          <xsl:for-each select="/nml:neuroml/net:populations/net:population[@name = $src]/net:instances/net:instance[@id = $preCellId]/net:location">
                              <xsl:variable name="zrotation">
                                <xsl:for-each select="../meta:notes">
                                  <xsl:if test="contains(text(),'zrotation')"> <!-- does text in the notes contain zrotation? -->
                                    <xsl:value-of select="substring-after(.,'zrotation=')"/>
                                  </xsl:if>
                                </xsl:for-each>
                              </xsl:variable>
                              <xsl:call-template name="writePoint">
                                <xsl:with-param name="CellType" select="$preCellType"/>
                                <xsl:with-param name="SegId" select="$preSegmentId"/>
                                <xsl:with-param name="Cellx" select="@x"/>
                                <xsl:with-param name="Celly" select="@y"/>
                                <xsl:with-param name="Cellz" select="@z"/>
                                <xsl:with-param name="zrotation" select="normalize-space($zrotation)"/>
                              </xsl:call-template>
                          </xsl:for-each>
                        </xsl:variable>
                        
                        <xsl:variable name="postPoint">
                          <xsl:for-each select="/nml:neuroml/net:populations/net:population[@name = $tgt]/net:instances/net:instance[@id = $postCellId]/net:location">
                              <xsl:variable name="zrotation">
                                <xsl:for-each select="../meta:notes">
                                  <xsl:if test="contains(text(),'zrotation')"> <!-- does text in the notes contain zrotation? -->
                                    <xsl:value-of select="substring-after(.,'zrotation=')"/>
                                  </xsl:if>
                                </xsl:for-each>
                              </xsl:variable>
                              <xsl:call-template name="writePoint">
                                <xsl:with-param name="CellType" select="$postCellType"/>
                                <xsl:with-param name="SegId" select="$postSegmentId"/>
                                <xsl:with-param name="Cellx" select="@x"/>
                                <xsl:with-param name="Celly" select="@y"/>
                                <xsl:with-param name="Cellz" select="@z"/>
                                <xsl:with-param name="zrotation" select="normalize-space($zrotation)"/>
                              </xsl:call-template>
                          </xsl:for-each>
                        </xsl:variable>
                        
                        <xsl:if test="($prePoint!='') and ($postPoint!='')">
                        <!-- draw all connections if drawGreyCells, else draw only two connections, color from setColorByConn below -->
                        <xsl:if test="$drawGreyCells or (($postCellType='granule') and (($postCellId=444) or ($postCellId=1160)))">

                          <xsl:variable name="colorRGBA">
                            <!--<xsl:call-template name="setColor">
                                <xsl:with-param name="colorIndex" select="$colorIndex"/>
                                <xsl:with-param name="subColorIndex" select="position()"/>                    
                            </xsl:call-template>-->
                            <xsl:call-template name="setColorByConn">
                                <xsl:with-param name="projName" select="$projName"/>
                                <xsl:with-param name="preCellType" select="$preCellType"/>
                                <xsl:with-param name="postCellType" select="$postCellType"/>
                                <xsl:with-param name="preCellId" select="$preCellId"/>
                                <xsl:with-param name="postCellId" select="$postCellId"/>                    
                            </xsl:call-template>
                          </xsl:variable>

                        <Transform>
                            <Shape>
                                <Appearance>
                                    <Material/>
                                        <xsl:if test="($preCellType='mitral') and ($postCellType='granule') and (($postCellId=444) or ($postCellId=1160))">
                                            <!-- Set the line width large for the colored connections. -->
                                            <xsl:if test="$setLineWidth='true'">
                                                <xsl:element name="LineProperties">
                                                    <xsl:attribute name="applied">true</xsl:attribute>
                                                    <xsl:attribute name="linetype">0</xsl:attribute>
                                                    <xsl:attribute name="linewidthScaleFactor"><xsl:value-of select="5"/></xsl:attribute>
                                                </xsl:element>
                                            </xsl:if>
                                        </xsl:if>
                                </Appearance>
                                <LineSet vertexCount="2">
                                    <xsl:element name="Coordinate">
                                        <xsl:attribute name="point"><xsl:value-of select="$prePoint"/>, <xsl:value-of select="$postPoint"/></xsl:attribute>
                                    </xsl:element>
                                    <!--<ColorRGBA color="1 0 0 0.3, 1 1 0 0.3"/>--><!-- red to yellow -->
                                    <!--<ColorRGBA color="0.5 0.5 0.5 0.3, 1 1 1 0.3"/>--><!-- grey to black -->
                                  <xsl:element name="ColorRGBA">
                                    <xsl:attribute name="color"><xsl:value-of select="$colorRGBA"/></xsl:attribute>
                                  </xsl:element>
                                </LineSet>
                            </Shape>
                        </Transform>
                        </xsl:if>
                        </xsl:if>

                    </xsl:if>
                </xsl:for-each> 
            </xsl:if>
        </xsl:for-each>

      </xsl:when>
      <!-- this displays single untranslated cells if there are no populations -->
      <xsl:otherwise>
        <xsl:for-each select="mml:cells/mml:cell | nml:cells/nml:cell">
            
            <xsl:variable name="cell_name"><xsl:value-of select="@name"/></xsl:variable>
            <xsl:call-template name="cell_draw">
                <xsl:with-param name="cell_type" select="$cell_name"/>
            </xsl:call-template>
            
        </xsl:for-each>
      </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template name="writePoint">
    <xsl:param name="CellType"/>
    <xsl:param name="SegId"/>
    <xsl:param name="Cellx"/>
    <xsl:param name="Celly"/>
    <xsl:param name="Cellz"/>
    <xsl:param name="zrotation"/>
    <xsl:for-each select="/nml:neuroml/nml:cells/nml:cell[@name = $CellType]/mml:segments/mml:segment[@id = $SegId]">
      <xsl:variable name="SegPx" select="mml:proximal/@x"/>
      <xsl:variable name="SegPy" select="mml:proximal/@y"/>
      <xsl:variable name="SegPz" select="mml:proximal/@z"/>
      <xsl:variable name="SegDx" select="mml:distal/@x"/>
      <xsl:variable name="SegDy" select="mml:distal/@y"/>
      <xsl:variable name="SegDz" select="mml:distal/@z"/>
      <xsl:variable name="Segx" select="($SegPx+$SegDx) div 2.0"/>
      <xsl:variable name="Segy" select="($SegPy+$SegDy) div 2.0"/>
      <xsl:variable name="Segz" select="($SegPz+$SegDz) div 2.0"/>
      <xsl:choose>
        <xsl:when test="$zrotation!=''"> <!-- if zrotation is provided -->
          <xsl:variable name="coszrot">
           <xsl:call-template name="cos">
              <xsl:with-param name="pX" select="$zrotation"/>
              <xsl:with-param name="pUnit" select="'rad'"/>
           </xsl:call-template>
          </xsl:variable>
          <xsl:variable name="sinzrot">
           <xsl:call-template name="sin">
              <xsl:with-param name="pX" select="$zrotation"/>
              <xsl:with-param name="pUnit" select="'rad'"/>
           </xsl:call-template>
          </xsl:variable>
          <xsl:variable name="RotatedSegx" select="$Segx*$coszrot - $Segy*$sinzrot"/>
          <xsl:variable name="RotatedSegy" select="$Segx*$sinzrot + $Segy*$coszrot"/>
          <xsl:variable name="RotatedSegz" select="$Segz"/>
          <!-- write the final translated and rotated point -->
          <xsl:value-of select="$Cellx+$RotatedSegx"/><xsl:text>  </xsl:text><xsl:value-of select="$Celly+$RotatedSegy"/><xsl:text>  </xsl:text><xsl:value-of select="$Cellz+$RotatedSegz"/>
        </xsl:when>
        <xsl:otherwise> <!-- if zrotation is not provided -->
          <!-- write the final translated point -->
          <xsl:value-of select="$Cellx+$Segx"/><xsl:text>  </xsl:text><xsl:value-of select="$Celly+$Segy"/><xsl:text>  </xsl:text><xsl:value-of select="$Cellz+$Segz"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
</xsl:template>

<xsl:template name="cell_draw">
    <xsl:param name="cell_type"/> <!-- empty cell_type -->
    <xsl:param name="x" select="'0'"/>
    <xsl:param name="y" select="'0'"/>
    <xsl:param name="z" select="'0'"/>
    <xsl:param name="zrotation" select="'0'"/>
    <xsl:param name="colorRGBA" select="'0 0 0 0, 0 0 0 0'"/>
    <xsl:for-each select="/nml:neuroml/nml:cells/nml:cell[@name=$cell_type]">
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
          <xsl:variable name="dia" select="mml:distal/@diameter"/>
          
          <xsl:element name="Transform">
            <xsl:attribute name="translation"><xsl:value-of select="$x"/><xsl:text> </xsl:text><xsl:value-of select="$y"/><xsl:text> </xsl:text><xsl:value-of select="$z"/></xsl:attribute>
            <xsl:if test="($zrotation!=0) and ($zrotation!='')">
                <xsl:attribute name="rotation"><xsl:text>0 0 1 </xsl:text><xsl:value-of select="$zrotation"/></xsl:attribute>
            </xsl:if>
              <Shape>
                  <Appearance>
                      <Material/>
                      <!-- Set the line width as segment diameter - but works only in Octaga; works in freewrl 1.22.13 too. -->
                      <xsl:if test="$setLineWidth='true'">
                          <xsl:element name="LineProperties">
                            <xsl:attribute name="applied">true</xsl:attribute>
                            <xsl:attribute name="linetype">0</xsl:attribute>
                            <xsl:attribute name="linewidthScaleFactor"><xsl:value-of select="$dia"/></xsl:attribute>
                          </xsl:element>
                      </xsl:if>
                  </Appearance>
                  <LineSet vertexCount="2">
                      <xsl:element name="Coordinate">
                          <xsl:attribute name="point"><xsl:value-of select="$proximal"/>, <xsl:value-of select="$distal"/></xsl:attribute>
                      </xsl:element>
                      <xsl:element name="ColorRGBA">
                        <xsl:attribute name="color"><xsl:value-of select="$colorRGBA"/></xsl:attribute>
                      </xsl:element>
                  </LineSet>
              </Shape>
          </xsl:element>
          
      </xsl:for-each>
    </xsl:for-each>
</xsl:template>

<xsl:template match="net:networkml">
        
    <xsl:for-each select="net:populations/net:population">
        
        <xsl:for-each select="net:instances/net:instance">
            <xsl:variable name="location"><xsl:value-of select="net:location/@x"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@y"/><xsl:text>  </xsl:text><xsl:value-of select="net:location/@z"/></xsl:variable>
            <xsl:element name="Transform">
                <xsl:attribute name="translation"><xsl:value-of select="$location"/></xsl:attribute>
        
        <Shape>
            <Appearance>
              <Material diffuseColor="0 1 0"/>
            </Appearance>
            <xsl:element name="Sphere">
                <xsl:attribute name="radius"><xsl:value-of select="$defaultCellRadius"/></xsl:attribute>
            </xsl:element>
        </Shape>   
            
            </xsl:element>
        </xsl:for-each>
    </xsl:for-each>
        
    <xsl:for-each select="net:projections/net:projection">
        <xsl:variable name="src"><xsl:value-of select="net:source"/><xsl:value-of select="@source"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
        <xsl:variable name="tgt"><xsl:value-of select="net:target"/><xsl:value-of select="@target"/></xsl:variable> <!-- Only one of attr or sub element should be present-->
        
        <xsl:if test="$src!='files'">

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
                            <!--<ColorRGBA color="0 1 0 0.3, 1 0 0 0.3"/>--><!-- Green to red-->
                            <ColorRGBA color="0.5 0.5 0.5 0.3, 1 1 1 0.3"/><!-- grey to black -->
                        </LineSet>
                    </Shape>
                </Transform>

            </xsl:for-each>

        </xsl:if>

    </xsl:for-each>

</xsl:template>

<xsl:template name="setColor">
    <xsl:param name="colorIndex" select="0"/>
    <xsl:param name="subColorIndex" select="0"/>
    <xsl:variable name="colormod" select="$colorIndex mod 6"/>
    <xsl:variable name="shade" select="(($subColorIndex mod 2)+1) div 2"/>
    <xsl:variable name="color">
        <xsl:choose> <!-- six colors: R, G, B, C, M, Y based on colorIndex; in two shades based on subColorIndex-->
            <xsl:when test="$colormod=0"><xsl:value-of select="$shade"/><xsl:text> </xsl:text>0<xsl:text> </xsl:text>0</xsl:when>
            <xsl:when test="$colormod=1">0<xsl:text> </xsl:text><xsl:value-of select="$shade"/><xsl:text> </xsl:text>0</xsl:when>
            <xsl:when test="$colormod=2">0<xsl:text> </xsl:text>0<xsl:text> </xsl:text><xsl:value-of select="$shade"/></xsl:when>
            <xsl:when test="$colormod=3"><xsl:value-of select="1-$shade"/><xsl:text> </xsl:text>1<xsl:text> </xsl:text>1</xsl:when>
            <xsl:when test="$colormod=4">1<xsl:text> </xsl:text><xsl:value-of select="1-$shade"/><xsl:text> </xsl:text>1</xsl:when>
            <xsl:when test="$colormod=5">1<xsl:text> </xsl:text>1<xsl:text> </xsl:text><xsl:value-of select="1-$shade"/></xsl:when>
        </xsl:choose>
    </xsl:variable>
    <xsl:value-of select="$color"/><xsl:text>, </xsl:text><xsl:value-of select="$color"/>
</xsl:template>

<xsl:template name="setColorByName">
    <xsl:param name="popName" select="mitrals"/>
    <xsl:param name="cellId" select="0"/>
    <xsl:variable name="glomNum" select="floor($cellId div 2) mod 3"/>
    <xsl:variable name="sisNumShade" select="($cellId mod 2)*0.5"/>
    <xsl:variable name="colorRGBA">
        <xsl:choose> <!-- R, G, B for gloms 0 to 3 (and cycled), in two shades for two sisters each. Transparent mild gray for others -->
            <xsl:when test="$popName='mitrals'">
                <xsl:choose>
                    <!--<xsl:when test="$glomNum=0">1<xsl:text> </xsl:text><xsl:value-of select="$sisNumShade"/><xsl:text> </xsl:text><xsl:value-of select="$sisNumShade"/><xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$glomNum=1"><xsl:value-of select="$sisNumShade"/><xsl:text> </xsl:text>1<xsl:text> </xsl:text><xsl:value-of select="$sisNumShade"/><xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$glomNum=2"><xsl:value-of select="$sisNumShade"/><xsl:text> </xsl:text><xsl:value-of select="$sisNumShade"/><xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>-->
                    <xsl:when test="$cellId=0">1<xsl:text> </xsl:text>0<xsl:text> </xsl:text>0<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$cellId=1">1<xsl:text> </xsl:text>0<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$cellId=2">0<xsl:text> </xsl:text>1<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$cellId=3">0<xsl:text> </xsl:text>0<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$cellId=4">0<xsl:text> </xsl:text>1<xsl:text> </xsl:text>0<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$cellId=5">1<xsl:text> </xsl:text>1<xsl:text> </xsl:text>0<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:otherwise>0.5<xsl:text> </xsl:text>0.5<xsl:text> </xsl:text>0.5<xsl:text> </xsl:text>0.2</xsl:otherwise>
                </xsl:choose>
            </xsl:when>
            <xsl:when test="($popName='granules_joints') and ($cellId=1160)"> <!-- mit 2 to mit 0 -->
                0<xsl:text> </xsl:text>1<xsl:text> </xsl:text>0<xsl:text> 1</xsl:text>
            </xsl:when>
            <xsl:when test="($popName='granules_joints') and ($cellId=444)"> <!-- mit 4 to mit 1 -->
                0<xsl:text> </xsl:text>0<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text>
            </xsl:when>
            <xsl:when test="($popName='granules_joints')"> <!-- all other joints -->
                <xsl:text>0 1 1 0.6</xsl:text>
            </xsl:when>
            <xsl:when test="($popName='granules_multis')"> <!-- all multis -->
                <xsl:text>0 0 1 0.6</xsl:text>
            </xsl:when>
            <xsl:when test="($popName='granules_singles')"> <!-- all singles -->
                <xsl:text>1 0 1 0.6</xsl:text>
            </xsl:when>
            <xsl:when test="($popName='PGs')"> <!-- all PGs -->
                <xsl:text>1 0 1 0.6</xsl:text>
            </xsl:when>
            <xsl:otherwise>0.5<xsl:text> </xsl:text>0.5<xsl:text> </xsl:text>0.5<xsl:text> </xsl:text>0.6</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>
    <xsl:value-of select="$colorRGBA"/><xsl:text>, </xsl:text><xsl:value-of select="$colorRGBA"/>
</xsl:template>

<xsl:template name="setColorByConn">
    <xsl:param name="projName" select="None"/>
    <xsl:param name="preCellType" select="mitral"/>
    <xsl:param name="postCellType" select="granule"/>
    <xsl:param name="preCellId" select="0"/>
    <xsl:param name="postCellId" select="0"/>
    <xsl:variable name="postCellIdInt" select="number($postCellId)"/>                  
    <xsl:variable name="colorRGBA">
        <xsl:choose>
            <!-- granule_joint #1320-1330 are between mits 1 and 4; #1331-1340 are between mits 0 and 2. Color only for these. -->
            <xsl:when test="($preCellType='mitral') and ($postCellType='granule') and (($postCellIdInt=444) or ($postCellIdInt=1160))">
                <xsl:choose> <!-- color as per other mitral color -->
                    <xsl:when test="$preCellId=0">0<xsl:text> </xsl:text>1<xsl:text> </xsl:text>0<xsl:text> </xsl:text>1</xsl:when>
                    <xsl:when test="$preCellId=1">0<xsl:text> </xsl:text>0<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>
                    <xsl:when test="$preCellId=2">0<xsl:text> </xsl:text>1<xsl:text> </xsl:text>0<xsl:text> 1</xsl:text></xsl:when>
                    <!--<xsl:when test="$preCellId=3">1<xsl:text> </xsl:text>1<xsl:text> </xsl:text>0<xsl:text> 1</xsl:text></xsl:when>-->
                    <xsl:when test="$preCellId=4">0<xsl:text> </xsl:text>0<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>
                    <!--<xsl:when test="$preCellId=5">0<xsl:text> </xsl:text>1<xsl:text> </xsl:text>1<xsl:text> 1</xsl:text></xsl:when>-->
                    <xsl:otherwise>0.4<xsl:text> </xsl:text>0.4<xsl:text> </xsl:text>0.4<xsl:text> </xsl:text>0.1</xsl:otherwise>
                </xsl:choose>
            </xsl:when>
            <xsl:when test="$projName='mitral_granule_main_exc_joints'">
                <xsl:text>0.4 0.4 0.4 0.4, 0 1 1 0.4</xsl:text>
            </xsl:when>
            <xsl:when test="$projName='mitral_granule_main_exc_multis'">
                <xsl:text>0.4 0.4 0.4 0.4, 0 0 1 0.4</xsl:text>
            </xsl:when>
            <xsl:when test="$projName='mitral_granule_main_exc_singles'">
                <xsl:text>0.4 0.4 0.4 0.4, 1 0 1 0.4</xsl:text>
            </xsl:when>
            <xsl:when test="$projName='mitral_PG'">
                <xsl:text>0.4 0.4 0.4 0.4, 1 0 1 0.4</xsl:text>
            </xsl:when>
            <xsl:otherwise>0.4<xsl:text> </xsl:text>0.4<xsl:text> </xsl:text>0.4<xsl:text> </xsl:text>0.4</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>
    <xsl:value-of select="$colorRGBA"/><xsl:text>, </xsl:text><xsl:value-of select="$colorRGBA"/>
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
