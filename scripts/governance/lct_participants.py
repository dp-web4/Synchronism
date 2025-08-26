#!/usr/bin/env python3
"""
LCT-based Participant Management for Whitepaper Governance
Implements Linked Context Tokens for participant identity and trust
"""

import json
import hashlib
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from enum import Enum
from dataclasses import dataclass, field, asdict

class ParticipantType(Enum):
    HUMAN = "human"
    AI_CLAUDE = "ai_claude"
    AI_GPT = "ai_gpt"
    AI_DEEPSEEK = "ai_deepseek"
    AI_OTHER = "ai_other"

class AccessMethod(Enum):
    API = "api"
    EMAIL = "email"
    GIT = "git"
    LOCAL = "local"

@dataclass
class ParticipantHistory:
    """Track participant's governance history"""
    proposals_created: int = 0
    reviews_submitted: int = 0
    proposals_accepted: int = 0
    proposals_rejected: int = 0
    arbiter_decisions: int = 0
    last_active: Optional[str] = None
    total_cycles_participated: int = 0
    holds_for_counter: List[str] = field(default_factory=list)  # Proposal IDs being held

@dataclass
class TrustScores:
    """Multi-dimensional trust scores (T3 tensor simplified)"""
    proposal_quality: float = 0.5  # Quality of proposals (0-1)
    review_quality: float = 0.5    # Quality of reviews (0-1)
    timeliness: float = 1.0        # Response within timeout (0-1)
    consistency: float = 0.5       # Consistency across contributions (0-1)
    collaboration: float = 0.5     # Willingness to iterate/improve (0-1)
    
    def aggregate(self) -> float:
        """Calculate aggregate trust score"""
        weights = {
            'proposal_quality': 0.3,
            'review_quality': 0.3,
            'timeliness': 0.2,
            'consistency': 0.1,
            'collaboration': 0.1
        }
        return sum(getattr(self, k) * v for k, v in weights.items())

@dataclass
class LCT:
    """Linked Context Token - Participant identity and state"""
    
    # Core identity
    id: str  # Unique identifier (hash of initial params)
    name: str
    type: ParticipantType
    created_at: str
    
    # Access configuration
    access_method: AccessMethod
    access_endpoint: str  # API URL, email address, or git remote
    access_credentials: Optional[Dict[str, str]] = None  # API key, git secret, etc.
    
    # Availability
    timeout_seconds: int = 300  # 5 minutes default
    last_heartbeat: Optional[str] = None
    is_available: bool = True
    
    # Accumulated state
    history: ParticipantHistory = field(default_factory=ParticipantHistory)
    trust_scores: TrustScores = field(default_factory=TrustScores)
    
    # Governance state
    current_cycle_proposed: bool = False
    current_cycle_reviewed: bool = False
    
    @classmethod
    def create(cls, 
               name: str,
               participant_type: ParticipantType,
               access_method: AccessMethod,
               access_endpoint: str,
               timeout_seconds: int = 300,
               access_credentials: Optional[Dict[str, str]] = None) -> 'LCT':
        """Create a new LCT for a participant"""
        
        # Generate unique ID from initial parameters
        id_source = f"{name}:{participant_type.value}:{access_endpoint}:{datetime.now().isoformat()}"
        lct_id = hashlib.sha256(id_source.encode()).hexdigest()[:16]
        
        return cls(
            id=lct_id,
            name=name,
            type=participant_type,
            created_at=datetime.now().isoformat(),
            access_method=access_method,
            access_endpoint=access_endpoint,
            access_credentials=access_credentials,
            timeout_seconds=timeout_seconds
        )
    
    def update_heartbeat(self):
        """Update last activity timestamp"""
        self.last_heartbeat = datetime.now().isoformat()
        self.is_available = True
    
    def check_availability(self) -> bool:
        """Check if participant is still available based on timeout"""
        # Humans are assumed available unless explicitly marked unavailable
        if self.type == ParticipantType.HUMAN:
            # For now, assume human is not available for automated cycles
            # This can be changed when human interaction is implemented
            return False
        
        if not self.last_heartbeat:
            return True  # Never pinged, assume available
        
        last_time = datetime.fromisoformat(self.last_heartbeat)
        timeout_delta = timedelta(seconds=self.timeout_seconds)
        
        if datetime.now() - last_time > timeout_delta:
            self.is_available = False
            return False
        
        return True
    
    def record_proposal(self, accepted: bool = False):
        """Record that participant created a proposal"""
        self.history.proposals_created += 1
        if accepted:
            self.history.proposals_accepted += 1
        else:
            self.history.proposals_rejected += 1
        self.history.last_active = datetime.now().isoformat()
        self.current_cycle_proposed = True
        
        # Update trust scores based on outcome
        if accepted:
            self.trust_scores.proposal_quality = min(1.0, self.trust_scores.proposal_quality + 0.1)
        else:
            self.trust_scores.proposal_quality = max(0.0, self.trust_scores.proposal_quality - 0.05)
    
    def record_review(self, was_helpful: bool = True):
        """Record that participant submitted a review"""
        self.history.reviews_submitted += 1
        self.history.last_active = datetime.now().isoformat()
        self.current_cycle_reviewed = True
        
        # Update trust scores
        if was_helpful:
            self.trust_scores.review_quality = min(1.0, self.trust_scores.review_quality + 0.05)
        self.trust_scores.collaboration = min(1.0, self.trust_scores.collaboration + 0.02)
    
    def record_arbiter_decision(self):
        """Record that participant made an arbiter decision"""
        self.history.arbiter_decisions += 1
        self.history.last_active = datetime.now().isoformat()
        self.trust_scores.consistency = min(1.0, self.trust_scores.consistency + 0.05)
    
    def add_hold_for_counter(self, proposal_id: str):
        """Mark that this participant wants to counter-propose"""
        if proposal_id not in self.history.holds_for_counter:
            self.history.holds_for_counter.append(proposal_id)
    
    def remove_hold(self, proposal_id: str):
        """Remove a hold after counter-proposal submitted"""
        if proposal_id in self.history.holds_for_counter:
            self.history.holds_for_counter.remove(proposal_id)
    
    def reset_cycle_flags(self):
        """Reset per-cycle activity flags"""
        self.current_cycle_proposed = False
        self.current_cycle_reviewed = False
        self.history.total_cycles_participated += 1
    
    def to_dict(self) -> Dict:
        """Convert LCT to dictionary for serialization"""
        return {
            'id': self.id,
            'name': self.name,
            'type': self.type.value,
            'created_at': self.created_at,
            'access_method': self.access_method.value,
            'access_endpoint': self.access_endpoint,
            'access_credentials': self.access_credentials,
            'timeout_seconds': self.timeout_seconds,
            'last_heartbeat': self.last_heartbeat,
            'is_available': self.is_available,
            'history': asdict(self.history),
            'trust_scores': asdict(self.trust_scores),
            'current_cycle_proposed': self.current_cycle_proposed,
            'current_cycle_reviewed': self.current_cycle_reviewed
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'LCT':
        """Restore LCT from dictionary"""
        lct = cls(
            id=data['id'],
            name=data['name'],
            type=ParticipantType(data['type']),
            created_at=data['created_at'],
            access_method=AccessMethod(data['access_method']),
            access_endpoint=data['access_endpoint'],
            access_credentials=data.get('access_credentials'),
            timeout_seconds=data.get('timeout_seconds', 300),
            last_heartbeat=data.get('last_heartbeat'),
            is_available=data.get('is_available', True)
        )
        
        # Restore history
        if 'history' in data:
            lct.history = ParticipantHistory(**data['history'])
        
        # Restore trust scores
        if 'trust_scores' in data:
            lct.trust_scores = TrustScores(**data['trust_scores'])
        
        lct.current_cycle_proposed = data.get('current_cycle_proposed', False)
        lct.current_cycle_reviewed = data.get('current_cycle_reviewed', False)
        
        return lct

class LCTRegistry:
    """Manages all participant LCTs"""
    
    def __init__(self, base_path: Optional[str] = None):
        if base_path is None:
            # Use relative path from current script location
            import os
            base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        self.base_path = Path(base_path)
        self.registry_file = self.base_path / "scripts" / "governance" / "config" / "lct_registry.json"
        self.participants: Dict[str, LCT] = {}
        self.load_registry()
    
    def load_registry(self):
        """Load LCTs from persistent storage"""
        if self.registry_file.exists():
            with open(self.registry_file, 'r') as f:
                data = json.load(f)
                for lct_id, lct_data in data.get('participants', {}).items():
                    self.participants[lct_id] = LCT.from_dict(lct_data)
    
    def save_registry(self):
        """Save LCTs to persistent storage"""
        self.registry_file.parent.mkdir(parents=True, exist_ok=True)
        
        data = {
            'participants': {
                lct_id: lct.to_dict() 
                for lct_id, lct in self.participants.items()
            },
            'last_updated': datetime.now().isoformat()
        }
        
        with open(self.registry_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def register_participant(self, lct: LCT) -> str:
        """Register a new participant"""
        self.participants[lct.id] = lct
        self.save_registry()
        return lct.id
    
    def get_participant(self, lct_id: str) -> Optional[LCT]:
        """Get participant by LCT ID"""
        return self.participants.get(lct_id)
    
    def get_available_participants(self, participant_type: Optional[ParticipantType] = None) -> List[LCT]:
        """Get all available participants, optionally filtered by type"""
        available = []
        for lct in self.participants.values():
            if lct.check_availability():
                if participant_type is None or lct.type == participant_type:
                    available.append(lct)
        return available
    
    def get_available_arbiters(self) -> List[LCT]:
        """Get available participants who can serve as arbiters"""
        arbiters = []
        for lct in self.participants.values():
            # Special handling for human arbiters
            if lct.type == ParticipantType.HUMAN:
                # Check if human is explicitly available (would need manual flagging)
                if getattr(lct, 'manually_available', False):
                    arbiters.append(lct)
                continue
            
            # For AI participants, check availability
            if not lct.check_availability():
                continue
            
            # High trust score can be arbiter
            if lct.trust_scores.aggregate() > 0.7:
                arbiters.append(lct)
            
            # Participants with enough history
            elif lct.history.total_cycles_participated > 5:
                arbiters.append(lct)
            
            # For initial cycles, allow any available AI to be arbiter
            elif lct.history.total_cycles_participated == 0:
                arbiters.append(lct)
        
        return arbiters
    
    def select_arbiter_with_fallback(self, preferred_id: Optional[str] = None) -> Optional[LCT]:
        """Select arbiter with timeout fallback"""
        
        # Try preferred arbiter first
        if preferred_id:
            arbiter = self.get_participant(preferred_id)
            if arbiter:
                # Check if it's a human who might be available
                if arbiter.type == ParticipantType.HUMAN:
                    if getattr(arbiter, 'manually_available', False):
                        return arbiter
                # Check if it's an available AI
                elif arbiter.check_availability():
                    return arbiter
        
        # Fallback to available arbiters
        available = self.get_available_arbiters()
        if not available:
            # Last resort: any available AI participant
            print("Warning: No qualified arbiters, selecting any available AI")
            for lct in self.participants.values():
                if lct.type != ParticipantType.HUMAN and lct.check_availability():
                    return lct
            return None
        
        # Sort by preference
        def arbiter_sort_key(lct):
            # Prefer certain AI types for arbitration
            type_preference = {
                ParticipantType.AI_CLAUDE: 0,  # Claude preferred for nuanced decisions
                ParticipantType.AI_GPT: 1,
                ParticipantType.AI_DEEPSEEK: 2,
                ParticipantType.AI_OTHER: 3,
                ParticipantType.HUMAN: 4  # Human last (usually unavailable)
            }
            
            return (
                type_preference.get(lct.type, 5),  # Type preference
                -lct.trust_scores.aggregate(),     # Higher trust better
                lct.history.arbiter_decisions,     # Fewer decisions (spread load)
                lct.history.last_active or "0"     # Least recently active
            )
        
        available.sort(key=arbiter_sort_key)
        
        selected = available[0]
        print(f"Selected arbiter: {selected.name} (trust: {selected.trust_scores.aggregate():.2f})")
        return selected
    
    def get_participants_for_cycle(self, max_participants: int = 10) -> List[LCT]:
        """Get participants eligible for current cycle"""
        eligible = []
        
        for lct in self.participants.values():
            if not lct.check_availability():
                continue
            
            # Skip if already participated this cycle
            if lct.current_cycle_proposed and lct.current_cycle_reviewed:
                continue
            
            eligible.append(lct)
            
            if len(eligible) >= max_participants:
                break
        
        return eligible
    
    def reset_cycle(self):
        """Reset all participants for new cycle"""
        for lct in self.participants.values():
            lct.reset_cycle_flags()
        self.save_registry()

def initialize_default_participants():
    """Initialize the registry with default participants"""
    registry = LCTRegistry()
    
    # Create LCT for human participant
    human_lct = LCT.create(
        name="Dennis",
        participant_type=ParticipantType.HUMAN,
        access_method=AccessMethod.EMAIL,
        access_endpoint="dennis@example.com",
        timeout_seconds=86400  # 24 hours for human
    )
    registry.register_participant(human_lct)
    
    # Create LCT for Claude
    claude_lct = LCT.create(
        name="Claude-3.5",
        participant_type=ParticipantType.AI_CLAUDE,
        access_method=AccessMethod.API,
        access_endpoint="https://api.anthropic.com/v1/messages",
        timeout_seconds=300,
        access_credentials={"api_key": "placeholder_key"}
    )
    registry.register_participant(claude_lct)
    
    # Create LCT for GPT
    gpt_lct = LCT.create(
        name="GPT-4",
        participant_type=ParticipantType.AI_GPT,
        access_method=AccessMethod.API,
        access_endpoint="https://api.openai.com/v1/chat/completions",
        timeout_seconds=300,
        access_credentials={"api_key": "placeholder_key"}
    )
    registry.register_participant(gpt_lct)
    
    # Create LCT for Deepseek
    deepseek_lct = LCT.create(
        name="Deepseek",
        participant_type=ParticipantType.AI_DEEPSEEK,
        access_method=AccessMethod.API,
        access_endpoint="https://api.deepseek.com/v1/chat",
        timeout_seconds=300,
        access_credentials={"api_key": "placeholder_key"}
    )
    registry.register_participant(deepseek_lct)
    
    print(f"Initialized {len(registry.participants)} participants:")
    for lct_id, lct in registry.participants.items():
        print(f"  - {lct.name} ({lct.type.value}): {lct_id}")
    
    return registry

if __name__ == "__main__":
    # Initialize default participants
    registry = initialize_default_participants()
    
    # Test availability and arbiter selection
    print("\nTesting arbiter selection with fallback:")
    arbiter = registry.select_arbiter_with_fallback("nonexistent_id")
    if arbiter:
        print(f"Selected arbiter: {arbiter.name} (trust: {arbiter.trust_scores.aggregate():.2f})")
    
    # Test cycle participants
    print("\nEligible participants for cycle:")
    participants = registry.get_participants_for_cycle(max_participants=3)
    for p in participants:
        print(f"  - {p.name}: proposed={p.current_cycle_proposed}, reviewed={p.current_cycle_reviewed}")